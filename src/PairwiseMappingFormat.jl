"""
    module PairwiseMappingFormat

Parse files in the Pairwise mApping Format (PAF).
"""
module PairwiseMappingFormat

export PAFReader,
    PAFRecord, aln_identity, aux_data, try_next!, query_coverage, target_coverage, is_mapped

# Not exported, because this very PAF-specific, and operates on foreign types.
public try_parse, ParserException, UnmappedError, Errors, Err

using PrecompileTools: @compile_workload, @setup_workload, @recompile_invalidations

# Note: I don't like to do this. Loading StringViews should cause widespread
# invalidation, but it does.
# This is because Base contain type unstable inference with AbstractString,
# which doesn't cause a lot of issues in Base since there aren't a lot of
# AbstractStrings in Base.
@recompile_invalidations begin
    using StringViews: StringView
end

using MemoryViews: ImmutableMemoryView

using XAMAuxData.SAM: Auxiliary
using BufferIO: AbstractBufReader, BufReader, ExpandReader, get_buffer, fill_buffer, consume

"""
    PAFRecord

Immutable type representing a single PAF line. The content of the record
is accessed through its properties.

A `PAFRecord` does not necessarily own its data, but may be a view into the
memory it was parsed from, e.g. the buffer of a [`PAFReader`](@ref) with `copy=false`.
Use `copy` to obtain a record that owns its data.

# Examples
```jldoctest
julia> record = PAFReader(first, open(path_to_paf));

julia> record.qname
"query1"

julia> record.qlen
301156
```

See also: [`PAFReader`](@ref)

# Extended help
The following properties may be used. For unmapped records, only `qname` and `qlen`
are available, and accessing any other property throws an [`UnmappedError`](@ref).
Use [`is_mapped`](@ref) to check if a record is mapped.
* `qname::StringView`. The query name, May be empty, and contain any bytes
  except `\\t` and `\\n`.
* `tname::StringView`. Target name. Like `qname`.
* `qlen` and `tlen::Int`, and gives the length of the query and target
  sequences, respectively. This must be > 0, and at most `typemax(UInt32)` (about 4.3 Gbp).
* `qstart`, `qend`, `tstart` and `tend::Int`, and give the starting and
  ending positions of the alignments on the query and target strand.
  These uses one-based, closed interval indexing as usual in Julia, but unlike
  the underlying PAF format.
  The intervals may be empty, in which case the ending position is the starting
  position minus one, so `qstart:qend` is an empty range. Hence, the ending positions
  are always >= the starting ones minus one, and always <= the query/target lengths.
* `matches::Int` gives the number of nucleotides / residues which match
  (i.e. are equal) in the alignment.
* `alnlen::Int` gives the length of the alignment
* `mapq::Union{Int, Nothing}` gives the mapping quality, or `nothing` if
  this information is unavailable.
* `is_rc::Bool` tells if the query and target strands are reverse-complement
  relative to each other.
"""
struct PAFRecord
    # The line the record was parsed from, without the trailing newline.
    # May be a view into memory not owned by the record, e.g. a reader's buffer.
    data::ImmutableMemoryView{UInt8}
    qname_len::Int32 # query name begins at index 1
    tname_start::Int32
    tname_len::Int32 # always zero when unmapped
    aux_start::Int32 # lastindex(data) + 1 if there is no aux data
    # Positions and lengths are stored as UInt32, so sequences of up to
    # typemax(UInt32) bytes are supported.
    qlen::UInt32 # always >= 1
    # Note: qstart/qend (and tstart/tend) are zero-indexed semi-open intervals
    # in the PAF format, and are stored as such. The properties are one-based.
    # The start is stored zero-based, such that the conversion can't overflow.
    # The fields below, and tname, are arbitrary when unmapped, and not accessible
    # through the properties of the record.
    qstart::UInt32
    qend::UInt32 # always >= qstart
    tlen::UInt32 # always >= 1
    tstart::UInt32
    tend::UInt32 # always >= tstart
    matches::UInt32
    alnlen::UInt32 # may be zero, and may be smaller than matches
    mapq::UInt8 # stored as 0xff for missing
    strand::UInt8 # 0x00: *, 0x01: +, 0x02: -
    # 16 bits of padding, which we can use for something else later
end

# Since the interface of the records are defined in terms of its
# properties, we need to be able to abstract over them
# Note: This is inlined so `sym` is constant folded, leaving only the branch for
# that property, and the cheap check for mapped-ness.
@inline function Base.getproperty(record::PAFRecord, sym::Symbol)
    return if sym === :qname
        qname(record)
    elseif sym === :qlen
        Int(getfield(record, :qlen))
    elseif sym === :tname
        check_mapped(record, sym)
        tname(record)
    elseif sym === :mapq
        check_mapped(record, sym)
        mapq(record)
    elseif sym in (:qstart, :tstart)
        check_mapped(record, sym)
        Int(getfield(record, sym)) + 1
    elseif sym in (:qend, :tlen, :tend, :matches, :alnlen)
        check_mapped(record, sym)
        Int(getfield(record, sym))
    elseif sym === :is_rc
        check_mapped(record, sym)
        is_rc(record)
    else
        error(lazy"Type PAFRecord has no property $sym")
    end
end

"""
    is_mapped(x::PAFRecord) -> Bool

Compute whether the `PAFRecord` is mapped.
Only the properties `qname` and `qlen`, and the auxiliary data of an unmapped record
are available. Accessing any other property of an unmapped record, or calling
[`query_coverage`](@ref), [`target_coverage`](@ref) or [`aln_identity`](@ref) on it,
throws an [`UnmappedError`](@ref).
Note that some PAF parsers do not handle unmapped records correctly, so be
wary when writing unmapped records.

# Examples
```jldoctest
julia> is_mapped(record)
true

julia> unmapped = parse(PAFRecord, "read1\\t1523\\t0\\t0\\t*\\t*\\t0\\t0\\t0\\t0\\t0\\t0");

julia> is_mapped(unmapped)
false
```

See also: [`PAFRecord`](@ref)
"""
is_mapped(x::PAFRecord) = !iszero(getfield(x, :strand))

"""
    UnmappedError

The exception thrown when accessing a property of an unmapped [`PAFRecord`](@ref) which
only exists for mapped records, such as `tname`, or when calling a function which
requires a mapped record, such as [`query_coverage`](@ref).
The name of the property or function may be accessed with the `.name` field.

# Examples
```jldoctest
julia> unmapped = parse(PAFRecord, "read1\\t1523\\t0\\t0\\t*\\t*\\t0\\t0\\t0\\t0\\t0\\t0");

julia> err = try unmapped.tname catch e e end;

julia> err isa PAF.UnmappedError
true

julia> err.name
:tname
```

See also: [`is_mapped`](@ref)
"""
struct UnmappedError <: Exception
    name::Symbol
end

function Base.showerror(io::IO, err::UnmappedError)
    return print(io, "UnmappedError: Cannot get `", err.name, "` of an unmapped PAFRecord. Check `is_mapped` first.")
end

# Keep the error path out of line, so the check is cheap to inline
@noinline throw_unmapped(name::Symbol) = throw(UnmappedError(name))

check_mapped(record::PAFRecord, name::Symbol) = (is_mapped(record) || throw_unmapped(name); nothing)

# Note: The functions below do not check if the record is mapped.
is_rc(record::PAFRecord) = getfield(record, :strand) == 0x02

# For tab completion
function Base.propertynames(::PAFRecord)
    return (
        :qname,
        :tname,
        :mapq,
        :qlen,
        :qstart,
        :qend,
        :tlen,
        :tstart,
        :tend,
        :matches,
        :alnlen,
        :is_rc,
    )
end

function PAFRecord()
    data = ImmutableMemoryView(Memory{UInt8}())
    return PAFRecord(data, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0xff, 0x00)
end

function Base.show(io::IO, ::MIME"text/plain", record::PAFRecord)
    buf = IOBuffer()
    is_mapped(record) || print(buf, "Unmapped ")
    println(buf, "PAFRecord:")
    println(buf, "  Query:    \"", qname(record), '"')
    if is_mapped(record)
        println(buf, "  Target:   \"", tname(record), '"')
        println(buf, "  Q cov:    ", round(query_coverage(record); digits = 4))
        println(buf, "  T cov:    ", round(target_coverage(record); digits = 4))
        println(buf, "  Identity: ", round(aln_identity(record); digits = 4))
        qual = mapq(record)
        qual = qual === nothing ? "255 (missing)" : string(qual)
        println(buf, "  Quality:  ", qual)
    end
    print(buf, "  Aux data: ")
    write(buf, repr_aux(record))
    return print(io, String(take!(buf)))
end

# Print the AUX fields indented
function repr_aux(record::PAFRecord)::Vector{UInt8}
    buf = IOBuffer()
    out = IOBuffer()
    show(buf, MIME"text/plain"(), aux_data(record))
    seekstart(buf)
    for line in eachline(buf; keep = true)
        write(out, "  ", line)
    end
    return take!(out)[3:end]
end

# Return StringView to not allocate
function qname(record::PAFRecord)
    return StringView(getfield(record, :data)[1:(getfield(record, :qname_len))])
end

# Return StringView to not allocate
function tname(record::PAFRecord)
    start = getfield(record, :tname_start)
    span = start:(start + getfield(record, :tname_len) - 1)
    return StringView(getfield(record, :data)[span])
end

# Note: We return nothing instead of missing because fuck missing
function mapq(record::PAFRecord)::Union{Int, Nothing}
    x = getfield(record, :mapq)
    return x == 0xff ? nothing : Int(x)
end

"""
    aln_identity(rec::Record) -> Float64

Return the nucleotide / residue identity, defined as the number of
matches divided by the alignment length.
For well-formed records, this is a number between 0 and 1. However, the parser does not
check that the number of matches is at most the alignment length, so this may be larger
than 1. If the alignment length is zero, this returns `NaN` or `Inf`.
Throws an [`UnmappedError`](@ref) if the record is unmapped.

# Examples
```jldoctest
julia> aln_identity(record)
1.0
```

See also: [`query_coverage`](@ref), [`target_coverage`](@ref)
"""
function aln_identity(record::PAFRecord)::Float64
    check_mapped(record, :aln_identity)
    return record.matches / record.alnlen
end

"""
    query_coverage(rec::Record) -> Float64

Compute the fraction of the query covered by the alignment.
Throws an [`UnmappedError`](@ref) if the record is unmapped.

# Examples
```jldoctest
julia> query_coverage(record)
0.9999535124653004
```

See also: [`target_coverage`](@ref), [`aln_identity`](@ref)
"""
function query_coverage(record::PAFRecord)::Float64
    check_mapped(record, :query_coverage)
    return (record.qend - record.qstart + 1) / record.qlen
end

"""
    target_coverage(rec::Record) -> Float64

Compute the fraction of the target covered by the alignment.
Throws an [`UnmappedError`](@ref) if the record is unmapped.

# Examples
```jldoctest
julia> target_coverage(record)
0.044934629307437725
```

See also: [`query_coverage`](@ref), [`aln_identity`](@ref)
"""
function target_coverage(record::PAFRecord)::Float64
    check_mapped(record, :target_coverage)
    return (record.tend - record.tstart + 1) / record.tlen
end

"""
    aux_data(rec::Record) -> XAMAuxData.SAM.Auxiliary

Return a lazily evaluated and validated `SAM.Auxiliary`, which is an `AbstractDict`
of the auxiliary fields at the end of the record.
For more details about this object, read the documentation of the package
XAMAuxData.jl.

# Examples
```jldoctest
julia> aux = aux_data(record);

julia> length(aux)
4

julia> aux["cm"]
29990

julia> haskey(aux, "k1")
false
```
"""
function aux_data(record::PAFRecord)
    return Auxiliary(getfield(record, :data)[getfield(record, :aux_start):end])
end

"""
    module Errors

Wrapper module used to contain the instances of error types,
and not pollute the namespace of the package.

The error kinds are intended for diagnostics. Which kind is reported for a given invalid
input may change in any minor release, and new kinds may be added.
See [`PairwiseMappingFormat.ParserException`](@ref).

# Examples
```jldoctest
julia> print(PAF.Errors.InvalidZero)
InvalidZero
```

See also: [`PairwiseMappingFormat.ParserException`](@ref)
"""
module Errors
    public Err
    # See descriptions in the `showerror` method below.
    @enum Err::Int32 begin
        TooFewFields
        IntegerOverflow
        BadInteger
        InvalidStrand
        EmptyNonTailingLine
        BackwardsIndices
        EmptyInteger
        InvalidZero
        PositionOutOfBounds
        LineTooLong
        EmbeddedNewline
    end
end # module Errors

using .Errors: Err

"""
    ParserException

The exception type thrown by iterating [`PAFReader`](@ref), or a failed `parse` of
a [`PAFRecord`](@ref). The functions `try_parse` and `try_next!` may return
(not throw) values of this type.

The line the error occurred on may be accessed with the `.line` property.
For `try_parse` and `parse`, this is always 1.

The `.kind` property is a [`PairwiseMappingFormat.Errors.Err`](@ref PairwiseMappingFormat.Errors)
describing why the line was rejected. It is intended for diagnostics, such as logging or
counting errors: An invalid line is often invalid in several ways, and which of them is
reported depends on implementation details. Therefore, the kind reported for a given
invalid input is deterministic within a version of this package, but may change in any
minor release, and new kinds may be added.
In contrast, whether input is rejected, and the line the error is reported on, are stable.

The error message printed by `showerror` includes the approximate position of the error
within the line. This position is not part of the API.

# Examples
```jldoctest
julia> err = PAF.try_parse("abc");

julia> err.line
1

julia> print(err.kind)
TooFewFields
```
"""
struct ParserException <: Exception
    kind::Errors.Err
    index_in_line::Int32
    line::Int
end

ParserException(index::Integer, kind::Errors.Err) = ParserException(kind, index, 1)

Base.propertynames(x::ParserException) = (:kind, :line)

function Base.getproperty(x::ParserException, s::Symbol)
    return if s ∈ (:kind, :line)
        getfield(x, s)
    else
        throw(ArgumentError(lazy"No such property in ParserException: $s"))
    end
end

function Base.showerror(io::IO, err::ParserException)
    buf = IOBuffer()
    print(buf, "Error when parsing PAF record on line ")
    print(buf, err.line)
    print(buf, ", near byte number ")
    print(buf, getfield(err, :index_in_line))
    print(buf, " in line: ")
    kind = err.kind
    s = if kind == Errors.TooFewFields
        "Not enough tab-separated fields. Each line must have at least 12 fields"
    elseif kind == Errors.IntegerOverflow
        "Overflow when parsing integer"
    elseif kind == Errors.BadInteger
        "Integer field contains characters out of 0-9 range"
    elseif kind == Errors.InvalidStrand
        "The strand column must be '+' or '-'"
    elseif kind == Errors.EmptyNonTailingLine
        "Whitespace-only line(s) followed by non-whitespace"
    elseif kind == Errors.BackwardsIndices
        "End coordinates must be equal to or larger than start coordinates for query and target"
    elseif kind == Errors.EmptyInteger
        "Integer field is empty"
    elseif kind == Errors.InvalidZero
        "Parsed an integer to zero which cannot meaningfully be zero"
    elseif kind == Errors.PositionOutOfBounds
        "Query or target start or end position is longer than query/target length"
    elseif kind == Errors.LineTooLong
        "Line is too long. PairwiseMappingFormat.jl can't handle lines of 2147483647 bytes or longer"
    elseif kind == Errors.EmbeddedNewline
        "Line contains a newline. Only a single trailing \"\\n\" or \"\\r\\n\" is allowed"
    end
    print(buf, s)
    return print(io, String(take!(buf)))
end

# Use this macro to return early and reduce code duplication.
# It's similar to `a = @something b return nothing`
macro var"?"(expr)
    return quote
        local res = $(esc(expr))
        res isa ParserException ? (return res) : res
    end
end

"""
    try_parse(x) -> Union{PAFRecord, ParserException}

Try parsing `x` into a `PAFRecord`. `x` may be any type that implements
`MemoryView`, such as `String` or `Vector{UInt8}`.

`x` must contain a single line, optionally ending with `"\\n"` or `"\\r\\n"`.
Any other `\\n` in `x` gives an error of kind `EmbeddedNewline`.
Other trailing whitespace, including a lone `\\r`, is not stripped:
It is part of the last field, and may make it invalid.

!!! warning
    The resulting record is a view into the memory of `x`, and does not copy it.
    If `x` is mutated, the record is mutated too, and may become invalid.
    Use `copy` on the record to obtain a record that owns its data.

# Examples
```jldoctest
julia> valid_line = open(readline, path_to_paf);

julia> typeof(PAF.try_parse(valid_line))
PAFRecord

julia> typeof(PAF.try_parse("invalid string"))
PairwiseMappingFormat.ParserException
```

See also: [`PAFRecord`](@ref), [`try_next!`](@ref), [`parse`](@ref Base.parse(::Type{PAFRecord}, ::Any))
"""
try_parse(x) = try_parse(ImmutableMemoryView(x))

# Strip one optional trailing \n or \r\n
function strip_newline(x::ImmutableMemoryView{UInt8})::ImmutableMemoryView{UInt8}
    (isempty(x) || @inbounds(x[end]) != UInt8('\n')) && return x
    has_cr = length(x) > 1 && @inbounds(x[end - 1]) == UInt8('\r')
    return @inbounds x[1:(end - 1 - has_cr)]
end

function try_parse(x::ImmutableMemoryView{UInt8})
    x = strip_newline(x)
    newline = findfirst(==(UInt8('\n')), x)
    newline === nothing || return ParserException(newline, Errors.EmbeddedNewline)
    return parse_line(x)
end

"""
    parse(::Type{PAFRecord}, x)::PAFRecord

Parse `x` into a `PAFRecord`, or throw a [`ParserException`](@ref
PairwiseMappingFormat.ParserException) if `x` is not a valid PAF line.
`x` may be any type that implements `MemoryView`, such as `String` or `Vector{UInt8}`.
See [`try_parse`](@ref PairwiseMappingFormat.try_parse) for what counts as a line.

!!! warning
    The resulting record is a view into the memory of `x`, and does not copy it.
    If `x` is mutated, the record is mutated too, and may become invalid.
    Use `copy` on the record to obtain a record that owns its data.

# Examples
```jldoctest
julia> record = parse(PAFRecord, codeunits("q\\t10\\t0\\t5\\t+\\tt\\t20\\t5\\t10\\t5\\t5\\t60\\n"));

julia> (record.qname, record.tname, record.mapq)
("q", "t", 60)
```
"""
function Base.parse(::Type{PAFRecord}, x)
    y = try_parse(x)
    return y isa PAFRecord ? y : throw(y)
end

function parse_line(mem::ImmutableMemoryView{UInt8})::Union{PAFRecord, ParserException}
    # Note: The aux data may begin at lastindex(mem) + 1, which must fit in an Int32
    lastindex(mem) >= typemax(Int32) &&
        return ParserException(typemax(Int32), Errors.LineTooLong)
    # Read fields incrementally. `i` is the index to read from, beginning at 1.
    # The parse_x functions will return the next index to begin from (after the tab)
    (qname, i) = @? parse_str(mem, 1)
    (qlen, i) = @? parse_int(mem, i, false)
    # Intervals may be empty, i.e. start == end, so the end may be zero
    (qstart, i) = @? parse_int(mem, i, true)
    (qend, i) = @? parse_int(mem, i, true)
    qend ≥ qstart || return ParserException(i % Int32, Errors.BackwardsIndices)
    qend > qlen && return ParserException(i % Int32, Errors.PositionOutOfBounds)

    # Load the strand field. We need at least 14 more bytes to encode 7 more mandatory fields
    # plus 7 more tabs.
    i + 14 > lastindex(mem) &&
        return ParserException(lastindex(mem) % Int32, Errors.TooFewFields)
    b = @inbounds mem[i]
    strand = if b == UInt8('*')
        0x00
    elseif b == UInt8('-')
        0x02
    elseif b == UInt8('+')
        0x01
    else
        return ParserException(i % Int32, Errors.InvalidStrand)
    end
    @inbounds mem[i + 1] == UInt8('\t') ||
        return ParserException((i + 1) % Int32, Errors.TooFewFields)
    iszero(strand) && return finish_unmapped(mem, qname, qlen, i + 2)

    # Load rest of the fields
    (tname, i) = @? parse_str(mem, i + 2) # i + 2 because strand plus tab took two bytes
    (tlen, i) = @? parse_int(mem, i, false)
    (tstart, i) = @? parse_int(mem, i, true)
    (tend, i) = @? parse_int(mem, i, true)
    tend ≥ tstart || return ParserException(i % Int32, Errors.BackwardsIndices)
    tend > tlen && return ParserException(i % Int32, Errors.PositionOutOfBounds)
    (matches, i) = @? parse_int(mem, i, true)
    (alnlen, i) = @? parse_int(mem, i, true)
    (mapq, i) = @? parse_int(mem, i, true, true)

    # A missing mapq is encoded as 255 in the PAF format, we store it as such
    mapq = if mapq > 255
        return ParserException(i % Int32, Errors.IntegerOverflow)
    else
        mapq % UInt8
    end

    # The aux data begins at `i`.
    # The integers have been checked to fit in an UInt32 by `parse_int`.
    return PAFRecord(
        mem,
        length(qname) % Int32,
        first(tname) % Int32,
        length(tname) % Int32,
        i % Int32,
        qlen % UInt32,
        qstart % UInt32,
        qend % UInt32,
        tlen % UInt32,
        tstart % UInt32,
        tend % UInt32,
        matches % UInt32,
        alnlen % UInt32,
        mapq,
        strand,
    )
end

# Parse the remaining fields of an unmapped record.
# These are validated like for mapped records, except that the integers may be zero,
# as written by e.g. minimap2's `--paf-no-hit`. Their values are not stored.
# i is the byte index of the strand plus two
function finish_unmapped(
        mem::ImmutableMemoryView{UInt8},
        qname::UnitRange{Int},
        qlen::Int,
        i::Int,
    )::Union{PAFRecord, ParserException}
    (_, i) = @? parse_str(mem, i)
    (tlen, i) = @? parse_int(mem, i, true)
    (tstart, i) = @? parse_int(mem, i, true)
    (tend, i) = @? parse_int(mem, i, true)
    tend ≥ tstart || return ParserException(i % Int32, Errors.BackwardsIndices)
    tend > tlen && return ParserException(i % Int32, Errors.PositionOutOfBounds)
    (_, i) = @? parse_int(mem, i, true) # matches
    (_, i) = @? parse_int(mem, i, true) # alnlen
    (mapq, i) = @? parse_int(mem, i, true, true)
    mapq > 255 && return ParserException(i % Int32, Errors.IntegerOverflow)

    # The aux data begins at `i`.
    return PAFRecord(mem, length(qname) % Int32, 1, 0, i % Int32, qlen % UInt32, 0, 0, 0, 0, 0, 0, 0, 0xff, 0x00)
end

# Just get the index of the next \t
function parse_str(
        v::ImmutableMemoryView{UInt8},
        from::Int,
    )::Union{Tuple{UnitRange{Int}, Int}, ParserException}
    i = findnext(==(UInt8('\t')), v, from)
    return if isnothing(i)
        ParserException(lastindex(v) % Int32, Errors.TooFewFields)
    else
        (from:(i - 1), i + 1)
    end
end

# TODO: 40% of runtime is spent in this function. Microoptimize it
function parse_int(
        v::ImmutableMemoryView{UInt8},
        from::Int,
        allow_zero::Bool,
        at_end::Bool = false,
    )::Union{Tuple{Int, Int}, ParserException}
    # Accumulate in an Int, which can't overflow because we check against
    # typemax(UInt32) after every digit. The result therefore fits in an UInt32.
    n = 0
    i = from
    for outer i in from:lastindex(v)
        b = @inbounds v[i]
        # If we hit a tab, this field is over.
        if b == UInt8('\t')
            i == from && return ParserException(i % Int32, Errors.EmptyInteger)
            (!allow_zero & iszero(n)) &&
                return ParserException(i % Int32, Errors.InvalidZero)
            return (n, i + 1)
        end
        b ∈ 0x30:0x39 || return ParserException(i % Int32, Errors.BadInteger)
        n = 10 * n + (b - 0x30)
        n > typemax(UInt32) && return ParserException(i % Int32, Errors.IntegerOverflow)
    end
    # Note: It is not possible to get an InvalidZero error here, because
    # the last field (mapq) does allow zero, and for any other fields,
    # if they do not end with a \t and thus reach this line, a TooFewFields
    # error will be returned.
    # at_end is if this is the mapq field, which does not need to end with a tab
    return if !at_end
        ParserException(lastindex(v) % Int32, Errors.TooFewFields)
    elseif from > lastindex(v)
        ParserException(i % Int32, Errors.EmptyInteger)
    else
        (n, i + 1)
    end
end

"""
    PAFReader(io::Union{IO, AbstractBufReader}; copy::Bool=true)

Construct a `PAFReader`, an iterator of [`PAFRecord`](@ref) read from `io`.
Iterating this object returns a `PAFRecord` for each valid line of input,
and throws a [`ParserException`](@ref) when an invalid record is found.

If `copy` is false, the records returned by iterating the reader, or by [`try_next!`](@ref),
are views into the reader's buffer, and are only valid until the reader is mutated,
which happens on iterating the reader or, `try_next!`, `isempty`.
This removes a few allocations and copies per iteration, but records must not be used
after they are invalidated.

See also: [`PAFRecord`](@ref), [`try_next!`](@ref)

# Examples
```jldoctest
julia> reader = PAFReader(open(path_to_paf)); typeof(reader)
PAFReader{BufferIO.BufReader{IOStream}}

julia> typeof(first(reader))
PAFRecord

julia> PAFReader(open(path_to_paf)) do reader
           for record in reader
                println(record.qlen)
           end
       end
301156
299273
288659

julia> PAFReader(open(path_to_paf); copy=false) do reader
           sum(record -> record.qlen, reader)
       end
889088
```

# Extended help
After throwing a `ParserException`, the reader advances over the malformed record,
so iteration may be consumed after catching the exception.

For efficiency, the auxiliary fields are not validated until they are accessed.

Lines end with `"\\n"` or `"\\r\\n"`. Other whitespace at the end of a line is not stripped:
It is part of the last field, and may make it invalid.

"""
mutable struct PAFReader{I <: AbstractBufReader}
    const io::ExpandReader{I}
    # Number of bytes of the last parsed line, including its newline, which have not
    # yet been consumed from `io`. Consuming is deferred to the next call to `try_next!`,
    # such that the record, which is a view into the buffer of `io`, stays valid.
    to_consume::Int
    line::Int
    # Whether to iterate copies of the parsed records
    copy::Bool
    # Whether a record which is a view into the buffer of `io` has been returned to the
    # user, and is still valid. While this is set, `isdone` must not fill the buffer.
    # Only ever set when `copy` is false.
    view_live::Bool
end

function PAFReader(io::AbstractBufReader; copy::Bool = true)
    return PAFReader{typeof(io)}(ExpandReader(io), 0, 1, copy, false)
end

function PAFReader(io::IO; buf_size::Int = 2^16, copy::Bool = true)
    return PAFReader(BufReader(io, max(buf_size, 16)); copy)
end

Base.IteratorSize(::Type{<:PAFReader}) = Base.SizeUnknown()
Base.eltype(::Type{<:PAFReader}) = PAFRecord
Base.close(reader::PAFReader) = close(reader.io)

function PAFReader(f::Function, io::Union{IO, AbstractBufReader}; kwargs...)
    reader = PAFReader(io; kwargs...)
    return try
        f(reader)
    finally
        close(reader)
    end
end

function Base.copy(record::PAFRecord)
    return PAFRecord(
        copy(getfield(record, :data)),
        getfield(record, :qname_len),
        getfield(record, :tname_start),
        getfield(record, :tname_len),
        getfield(record, :aux_start),
        getfield(record, :qlen),
        getfield(record, :qstart),
        getfield(record, :qend),
        getfield(record, :tlen),
        getfield(record, :tstart),
        getfield(record, :tend),
        getfield(record, :matches),
        getfield(record, :alnlen),
        getfield(record, :mapq),
        getfield(record, :strand),
    )
end

function Base.iterate(reader::PAFReader, (::Nothing) = nothing)
    res = try_next!(reader)
    return if res isa ParserException
        throw(res)
    elseif res === nothing
        nothing
    else
        (res, nothing)
    end
end

"""
    try_next!(reader::PAFReader) -> Union{Nothing, PAFRecord, ParserException}

Try to read a record from the [`PAFReader`](@ref) and advance the reader. This may yield one of three options:
* If the reader is empty, return `nothing` and do not advance the reader
* If the next record is invalid, return (do not throw) a `ParserException` and
  advance the reader past the invalid line. The next call continues from the next line,
  so invalid lines may be skipped by ignoring the error.
* If the next record is valid, return it as a [`PAFRecord`](@ref) and advance the reader

If the reader was constructed with `copy=true` (the default), the returned record owns
its data. If `copy=false`, the record is a view into the reader's buffer, which is
invalidated when the reader is advanced, and may be invalidated by calling
[`isempty`](@ref Base.isempty(::PAFReader)) on the reader. See [`PAFReader`](@ref).

# Examples
```jldoctest
julia> reader = PAFReader(open(path_to_paf));

julia> try_next!(reader) isa PAFRecord
true

julia> [try_next!(reader) for i in 1:2] isa Vector{PAFRecord}
true

julia> typeof([try_next!(reader) for i in 1:100])
Vector{Nothing} (alias for Array{Nothing, 1})

julia> close(reader); reader = PAFReader(IOBuffer("bad data\\n" * open(readline, path_to_paf)));

julia> try_next!(reader) isa PAF.ParserException
true

julia> try_next!(reader) isa PAFRecord
true

julia> try_next!(reader) === nothing
true
```
"""
function try_next!(reader::PAFReader)::Union{Nothing, PAFRecord, ParserException}
    res = next_view!(reader)
    res isa PAFRecord || return res
    # A copy owns its data, so the view never reaches the user
    reader.copy && return copy(res)
    reader.view_live = true
    return res
end

# Like `try_next!`, but always returns a view, which is not yet exposed to the user
function next_view!(reader::PAFReader)::Union{Nothing, PAFRecord, ParserException}
    consume_previous!(reader)
    (buffer, last_index, line_len) = buffer_line(reader.io)
    line = @inbounds buffer[1:last_index]
    res = parse_line(line)
    # A whitespace-only line, which can never be a valid record, is the end of the
    # stream if the rest of it is whitespace. Checking only on errors keeps the
    # check off the hot path.
    if res isa ParserException && all(is_whitespace, line)
        tail = check_tail(reader, buffer)
        tail === nothing && return nothing
        # Skip all the whitespace-only lines, and report the error at the first of them
        (n_bytes, n_lines) = tail
        err = ParserException(Errors.EmptyNonTailingLine, Int32(1), reader.line)
        reader.to_consume = n_bytes
        reader.line += n_lines
        return err
    end
    res = if res isa ParserException
        ParserException(res.kind, getfield(res, :index_in_line), reader.line)
    else
        res
    end
    # Advance past the line, even if it is invalid, such that the next call
    # continues from the next line
    reader.to_consume = line_len
    reader.line += 1
    return res
end

# Consume the previous line. This invalidates the record if it is a view into the buffer.
function consume_previous!(reader::PAFReader)
    if !iszero(reader.to_consume)
        consume(reader.io, reader.to_consume)
        reader.to_consume = 0
    end
    reader.view_live = false
    return nothing
end

# Fill the buffer until it contains a full line, or EOF is reached.
# Return the buffer, the index of the last byte of the line before the \n (or \r\n),
# and the number of bytes in the line including the newline.
function buffer_line(io::ExpandReader)::Tuple{ImmutableMemoryView{UInt8}, Int, Int}
    buffer = get_buffer(io)
    newline = findfirst(==(0x0a), buffer)
    # Keep reading until we find a newline or reach EOF
    while newline === nothing
        # The number of bytes in buffer we have already memchrd searching for a newline
        n_scanned = length(buffer)
        # This signals EOF - we might still have data in the buffer, if the underlying stream
        # ends with a record that doesn't end in a newline
        iszero(fill_buffer(io)) && return (buffer, length(buffer), length(buffer))
        buffer = get_buffer(io)
        newline = findnext(==(0x0a), buffer, n_scanned + 1)
    end
    # Last index before \n, and before \r if the newline is \r\n as on Windows
    has_windows_newline = !isone(newline) && @inbounds(buffer[newline - 1]) == UInt8('\r')
    return (buffer, newline - (1 + has_windows_newline), newline)
end

is_whitespace(b::UInt8) = b ∈ (UInt8('\r'), UInt8('\n'), UInt8('\t'), UInt8(' '))

# Called when the next line is whitespace-only, with `buffer` beginning at that line.
# Return `nothing` if the rest of the stream is whitespace. Else, return the number of
# bytes and lines of the whitespace-only lines before the first line with non-whitespace.
# The buffer may end after the line, so keep filling it until EOF.
# This does not mutate `reader`, except by filling its buffer.
function check_tail(reader::PAFReader, buffer::ImmutableMemoryView{UInt8})::Union{Nothing, Tuple{Int, Int}}
    io = reader.io
    n_lines = 0
    # Index of the first byte after the last newline seen
    line_start = 1
    n_scanned = 0
    while true
        for i in (n_scanned + 1):length(buffer)
            b = @inbounds buffer[i]
            if b == UInt8('\n')
                n_lines += 1
                line_start = i + 1
            elseif !is_whitespace(b)
                return (line_start - 1, n_lines)
            end
        end
        n_scanned = length(buffer)
        iszero(fill_buffer(io)) && return nothing
        buffer = get_buffer(io)
    end
    return
end

# This does not advance the reader, and never invalidates a returned record.
# If that is not possible, because a returned record is a view into the buffer and the
# buffer must be filled to determine the answer, return `missing`.
function Base.isdone(reader::PAFReader, _...)::Union{Bool, Missing}
    rest = get_buffer(reader.io)[(reader.to_consume + 1):end]
    # The reader is done exactly when the rest of the stream is whitespace:
    # If the next line contains non-whitespace, `try_next!` returns a record or an error.
    # Else, it is the end of the stream only if the rest of it is whitespace.
    # This usually returns at the first byte.
    all(is_whitespace, rest) || return false
    # If no view is live, we can freely consume the previous line and fill the buffer.
    # Else, we can't tell if the stream continues with non-whitespace.
    return reader.view_live ? missing : isdone_filling!(reader)
end

function isdone_filling!(reader::PAFReader)::Bool
    consume_previous!(reader)
    return check_tail(reader, get_buffer(reader.io)) === nothing
end

"""
    isempty(reader::PAFReader)::Bool

Check whether `reader` has no more records, i.e. whether [`try_next!`](@ref) would return
`nothing`. This does not advance the reader, but if the reader was constructed with
`copy=false`, it may invalidate the previous record by mutating the reader's inner state.

In contrast, `Base.isdone(reader)`, never invalidates a record,
and returns `missing` if it can't determine whether the reader is done without invalidating
previous records.
"""
function Base.isempty(reader::PAFReader)
    done = Base.isdone(reader)
    done === missing || return done
    # Invalidate the live view, such that `isdone` may fill the buffer
    return isdone_filling!(reader)
end

@setup_workload begin
    path = joinpath(pkgdir(PairwiseMappingFormat), "test", "example.paf")
    @compile_workload begin
        record = PAFReader(open(path)) do reader
            (rec, _) = iterate(reader)
            _ = try_next!(reader)
        end
        line = first(eachline(path))
        rec = parse(PAFRecord, line)
        ncodeunits(rec.qname) + copy(rec).qlen + rec.mapq + ncodeunits(rec.tname)
    end
end

end # module PairwiseMappingFormat
