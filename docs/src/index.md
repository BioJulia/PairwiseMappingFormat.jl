```@meta
CurrentModule = PairwiseMappingFormat
DocTestSetup = quote
    using PairwiseMappingFormat
    using StringViews: StringView
    path_to_paf = joinpath(pkgdir(PairwiseMappingFormat), "test", "example.paf")
    record = PAFReader(first, open(path_to_paf))
    const PAF = PairwiseMappingFormat
    using PairwiseMappingFormat: try_parse
end
```

# PairwiseMappingFormat.jl
This package is for reading files of the PAF (Pairwise mApping Format) format,
which is a simple tab-delimited format used by e.g. minimap2 and strobealign
to store records representing pairwise alignment between biological sequences.

## Reader
The [`PAFReader`](@ref) type wraps an `IO`, and is an iterator of [`PAFRecord`](@ref) objects:

```jldoctest; output = false
reader = PAFReader(open(path_to_paf))
records = collect(reader)
close(reader)

@assert isempty(reader)
@assert typeof(records) == Vector{PAFRecord}

# output

```

Similar to the common `open(path) do io`-syntax in Julia, `PAFReader` takes an optional
`f::Function` as a first argument, in which case it will apply `f` to the returned
`PAFReader` object, and close the underlying io when `f` returns (or errors):

```jldoctest
julia> PAFReader(open(path_to_paf)) do reader
           for record in reader
               println(record.qlen)
           end
       end
301156
299273
288659
```

The [`PAFReader`](@ref) constructor takes an optional keyword `copy`, which defaults to `true`.
If it is `false`, the record will iterate the _same_ [`PAFRecord`](@ref) object, overwriting
it in-place.
This reduces allocations and gives a slight speed increase, but may result in bugs
if records, or references of records of previous iterations are stored and
unexpectedly overwritten:

```jldoctest
julia> records = PAFReader(collect, open(path_to_paf); copy=false);

julia> println(map(i -> i.qlen, records)) # NB: All the same record!
[288659, 288659, 288659]
```

At the moments, readers do not support seeking.

## Record
The mutable [`PAFRecord`](@ref) object represents a single line in a PAF file.
The individual columns of the PAF line is obtained by accessing properties of the records:

```jldoctest record
julia> record = PAFReader(first, open(path_to_paf));

julia> println(record.qname)
query1

julia> println(record.qlen)
301156

julia> println(record.mapq)
0
```

Note that `Base.getproperty` is overloaded for `PAFRecord`, so the properties
are public and stable across major versions, but may not reflect the actual
underlying fields as they are stored in the `PAFRecord`.

### Fields of `PAFRecord`
* `qname::StringView{A}`. `A` is an implementation detail.
  Names of the query sequences
  When the sequences are from a FASTA file, this is typically the _identifier_ of the
  FASTA records, i.e. the part of the header before the first space.
* `tname::Union{Nothing, StringView{A}}`. Same as `qname`, but for the target (i.e. subject)
  sequence. May be `nothing` if the record is unmapped.
* `qlen` and `tlen`. Of type `Int`. Length of the full query and target sequence.
* `qstart` and `tstart`. Of type `Int`. Leftmost (i.e. lowest) position of the alignment
  in the query and target, respectively. This does not take strandedness into account.
* `qend` and `tend`. Of type `Int`. Rightmost (i.e. highest) position of the alignment
  in the query and target, respectively. This does not take strandedness into account.
  As such, we always have `record.qend - record.qtart) >= 0`, and likewise for the target.
* `alnlen::Int`. Length of alignment. That includes gaps in query and target, and may therefore
  not match either `qend - qstart + 1`, or that of the target.
* `matches::Int`. Number of residue matches (not mismatches) in the alignment.
  `matches / alnlen` is the alignment identity and is between 0 and 1.
* `mapq::Union{Int, Nothing}`. Mapping quality, from 0:254, where 254 is the better quality.
  When calibrated, the probability the mapping is correct should be 10^(-mapq / 10).
  A value of `nothing` indicates the mapping quality is unavailable.
* `is_rc::Union{Bool, Nothing}` indicates whether the alignment matches the target on the forward (`false`) or reverse-complement `true` strand. A missing alignment is `nothing`.

The PAF format does not contain fields for alignment identity, query coverage
or target coverage.
However, these can be approximated with the functions [`aln_identity`](@ref),
[`query_coverage`](@ref) and [`query_coverage`](@ref).

Besides being iterated from [`PAFReader`](@ref), records can also be created by
parsing from a bytes-like object such as a `String`:

```jldoctest; output = false
data = "contig_11\t288659\t27\t288653\t+\tCP004047.1\t6701780\t5027225\t5315851\t288618\t288626\t0\ttp:A:P\tcm:i:28871\ts1:i:288618\ts2:i:288618\tdv:f:0.0000\trl:i:57"

record = parse(PAFRecord, data)
@assert record isa PAFRecord

# output

```

Rarely, PAF records are unmapped. An unmapped record contain the query id, the query length, and any [auxiliary fields](@ref aux).
The target name and the strand of an unmapped return is `nothing`.
You can check if a record is mapped or unmapped with the function `is_mapped`:

```jldoctest unmapped
julia> unmapped = parse(PAFRecord, "myquery\t5032251\t9\t11\t*\t*\t5\t2\t1\t7\t4\t0");

julia> @assert !is_mapped(unmapped)

julia> println(unmapped.qname)
myquery

julia> println(unmapped.qlen)
5032251

julia> @assert isnothing(unmapped.tname)

julia> # The strand is given by is_rc, and is nothing for unmapped records
       @assert isnothing(unmapped.is_rc)
```

All other properties of unmapped records may contain arbitrary data, and should not be relied on.
For this reason, the functions [`target_coverage`](@ref), [`query_coverage`](@ref) and [`aln_identity`](@ref)
may give nonsense answers for unmapped records:
```jldoctest unmapped
julia> # Note! These results are arbitrary and not guaranteed to be stable
       # across minor releases of this package
       println(target_coverage(unmapped))
Inf

julia> println(query_coverage(unmapped))
1.9871822768776836e-7

julia> println(aln_identity(unmapped))
NaN
```

## Low-level interface
Iterating `PAFReader`s, and the `parse` function will throw a `PairwiseMappingFormat.ParserException` if the data is invalid:
```jldoctest
julia> parse(PAFRecord, "not a PAF line")
ERROR: Error when parsing PAF record on line 1, near byte number 14 in line: Not enough tab-separated fields. Each line must have at least 12 fields
[...]
```

In some circuomstances, throwing exceptions may not be acceptable, and so PairwiseMappingFormat.jl
provide functionality for returning errors as values.
These values can then be checked to programmatically handle error states.

The public, unexported [`try_parse`](@ref) function can be used instead of `Base.parse`.
It either returns a valid `PAFRecord`, or else returns the [`ParserException`](@ref),
but does not throw the exception:

```jldoctest public
julia> const PAF = PairwiseMappingFormat;

julia> println(PAF.try_parse("not a PAF line"))
PairwiseMappingFormat.ParserException(PairwiseMappingFormat.Errors.TooFewFields, 14, 1)
```

Similarly, the next record of a `PAFReader` may be obtained with the unexported [`try_next!`](@ref)
function.
This function will either:
* Return `nothing` if the reader has no more records
* Return a `ParserException` if the next record is invalid
* Return a `PAFRecord` if there is another valid record.
  Only in this case it will advance the position of the reader.

In other words, if a call to [`try_next!`](@ref) returns `nothing` or a `ParserException`,
any future calls will return the same value forever:

```jldoctest public; output = false
reader = PAFReader(open(path_to_paf))

# Emits `PAFRecord` for valid records
records = [try_next!(reader) for i in 1:3]
@assert records isa Vector{PAFRecord}

# Indefinitely emits `nothing` once the reader is empty
nothings = [try_next!(reader) for i in 1:10]
@assert nothings isa Vector{Nothing}

close(reader)

reader = PAFReader(IOBuffer("not a PAF record"))
err = try_next!(reader)

@assert err isa PAF.ParserException
err2 = try_next!(reader)
@assert err == err2

# output

```

The [`ParserException`](@ref) type contains the error type as an `Enum`, and the line where
the exception occurred. These can be obtained with the `.kind` and `.line` properties.

```jldoctest public
julia> err.line
1

julia> err.kind
TooFewFields::Err = 0
```

The precise value of this `.kind` object for any given error condition is subject
to change across minor versions and new values may be introduced.
This is because the same error may be detected in multiple different ways.

## [Auxiliary fields](@id aux)
[`PAFRecord`](@ref)s may contain extra auxiliary fields at the end of the records,
similar to SAM records.
Any auxiliary data is stored in the `PAFRecord`, and lazily parsed and validated
as they are accessed. This means that a `PAFRecord` may contain invalid auxiliary
data.

The [`aux_data`](@ref) function constructs a lazy `SAM.Auxiliary` object from the package
XAMAuxData.jl.
The precise semantics of this `Auxiliary` type is a little complicated, and
can be found in the documentation of XAMAuxData.jl.
Here is a small example:

```jldoctest
julia> line = first(eachline(path_to_paf));

julia> aux_string = join(split(line, '\t')[13:end], '\t')
"tp:A:P\tcm:i:29990\tdv:f:0.0000\tkm:Z:some text here"

julia> record = parse(PAFRecord, line);

julia> aux = aux_data(record);

julia> isvalid(aux)
true

julia> aux["tp"]
'P': ASCII/Unicode U+0050 (category Lu: Letter, uppercase)

julia> aux["km"]
"some text here"

julia> haskey(aux, "cm")
true

julia> aux["AB"] = 5.5;

julia> haskey(aux, "AB")
true
```

## Misc info
* The `PAFReader` cannot handle trailing whitespace, except trailing newlines at
  the end of the file. This is because trailing whitespace may be significant
  in some records, e.g. if it ends with an auxiliary field of the element type
  `Z`.
  A PAF record with trailing whitespace is considered invalid.


## Reference
```@docs
PairwiseMappingFormat.PairwiseMappingFormat
PAFReader
PAFRecord
aux_data
aln_identity
query_coverage
target_coverage
try_next!
is_mapped
PairwiseMappingFormat.Errors
PairwiseMappingFormat.try_parse
PairwiseMappingFormat.ParserException
```
