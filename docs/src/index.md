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
which is a simple tab-delimited format used by e.g. minimap2 and strobealign.

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
PAFReader(open(path_to_paf)) do reader
    for record in reader
        println(record.qlen)
    end
end

# output
301156
299273
288659

```

The [`PAFReader`](@ref) constructor takes an optional keyword `copy`, which defaults to `true`.
If it is `false`, the record will iterate the _same_ [`PAFRecord`](@ref) object, overwriting
it in-place.
This reduces allocations and give a slight speed increase, but may result in bugs
if records, or references of records of previous iterations are stored:

```jldoctest
records = PAFReader(collect, open(path_to_paf); copy=false)
println(map(i -> i.qlen, records))

# output
[288659, 288659, 288659]

```

At the moments, readers do not support seeking.

## Record
The mutable [`PAFRecord`](@ref) object represents a single line in a PAF file.
The individual columns of the PAF line is obtained by accessing properties of the records:

```jldoctest record
record = PAFReader(first, open(path_to_paf))

println(record.qname)
println(record.qlen)
println(record.mapq)

# output
query1
301156
0
```

Note that `Base.getproperty` is overloaded for `PAFRecord`, so the properties
are public and stable across major versions, but may not reflect the actual
underlying fields as they are stored in the `PAFRecord`.

A description of the properties of `PAFRecord`s can be found in the docstring of
`PAFRecord`, but here is a list of them:

```jldoctest record; output = false
properties = [
    (:qname,   StringView),
    (:qlen,    Int),
    (:qstart,  Int),
    (:qend,    Int),
    (:tname,   StringView),
    (:tlen,    Int),
    (:tstart,  Int),
    (:tend,    Int),
    (:matches, Int),
    (:alnlen,  Int),
    (:mapq,    Union{Int, Nothing}),
    (:is_rc,   Bool),
]

@assert length(propertynames(record)) == length(properties)

for (property, type) in properties
    @assert getproperty(record, property) isa type
end

# output

```

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

## Low-level interface
Iterating `PAFReader`s, and the `parse` function will throw a `PairwiseMappingFormat.ParserException` if the data is invalid:
```jldoctest
parse(PAFRecord, "not a PAF line")

# output
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
const PAF = PairwiseMappingFormat

println(PAF.try_parse("not a PAF line"))

# output
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
println(err.line)
println(err.kind)

# output
1
TooFewFields
```

The precise value of this `.kind` object for any given error condition is subject
to change across minor versions and new values may be introduced.
This is because the same error may be detected in multiple different ways.

## Auxiliary fields
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
PAFReader
PAFRecord
aux_data
aln_identity
query_coverage
target_coverage
try_next!
PairwiseMappingFormat.try_parse
PairwiseMappingFormat.ParserException
```
