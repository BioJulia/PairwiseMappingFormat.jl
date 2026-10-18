"""
    module PairwiseMappingFormat

Parse files of the PairwisemAppingFormat (PAF) format.
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

Mutable type representing a single PAF line. The content of the record
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
  sequences, respectively. This must be > 0.
* `qstart`, `qend`, `tstart` and `tend::Int`, and give the starting and
  ending positions of the alignments on the query and target strand.
  These uses one-based, closed interval indexing as usual in Julia, but unlike
  the underlying PAF format.
  The ending positions are always >= the starting ones, and the ending positions
  are always <= the query/target lengths.
* `matches::Int` gives the number of nucleotides / residues which match
  (i.e. are equal) in the alignment.
* `alnlen::Int` gives the length of the alignment
* `mapq::Union{Int, Nothing}` gives the mapping quality, or `nothing` if
  this information is unavailable.
* `is_rc::Bool` tells if the query and target strands are reverse-complement
  relative to each other.
"""
mutable struct PAFRecord
    # The line the record was parsed from, without the trailing newline.
    # May be a view into memory not owned by the record, e.g. a reader's buffer.
    data::ImmutableMemoryView{UInt8}
    qname_len::Int32 # query name begins at index 1
    tname_start::Int32
    tname_len::Int32 # always zero when unmapped
    aux_start::Int32 # lastindex(data) + 1 if there is no aux data
    qlen::Int32 # always >= 1
    # Note: qstart/qend (and tstart/tend) are zero-indexed semi-open interval
    # in PAF format, but stored in this struct one-based.
    # The fields below, and tname, are arbitrary when unmapped, and not accessible
    # through the properties of the record.
    qstart::Int32 # always >= 1
    qend::Int32 # always >= qstart
    tlen::Int32 # always >= 1
    tstart::Int32 # always >= 1
    tend::Int32 # always >= tstart
    matches::Int32
    alnlen::Int32 # always >= 1
    mapq::UInt8 # stored as 0xff for missing
    strand::UInt8 # 0x00: *, 0x01: -, 0x02: +
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
    elseif sym in (:qstart, :qend, :tlen, :tstart, :tend, :matches, :alnlen)
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
    return PAFRecord(data, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0xff, 0x00)
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
This is a number between 0 and 1.
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
    return (record.qend - record.qstart + Int32(1)) / record.qlen
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
    return (record.tend - record.tstart + Int32(1)) / record.tlen
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
    end
end # module Errors

using .Errors: Err

"""
    ParserException

The exception type thrown by iterating [`PAFReader`](@ref), or a failed `parse` of
a [`PAFRecord`](@ref). The functions `try_parse` and `try_next!` may return
(not throw) values of this type.
The line the error occurred may be accessed with the `.line` field.
The error kind may be accessed with the `.kind` field.

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

See also: [`PAFRecord`](@ref), [`try_next!`](@ref)
"""
try_parse(x) = try_parse(ImmutableMemoryView(x))
try_parse(x::ImmutableMemoryView{UInt8}) = parse_line!(PAFRecord(), x)

function Base.parse(::Type{PAFRecord}, s::Union{String, SubString{String}})
    y = try_parse(s)
    return y isa PAFRecord ? y : throw(y)
end

function parse_line!(
        record::PAFRecord,
        mem::ImmutableMemoryView{UInt8},
    )::Union{PAFRecord, ParserException}
    # Note: The aux data may begin at lastindex(mem) + 1, which must fit in an Int32
    lastindex(mem) >= typemax(Int32) &&
        return ParserException(typemax(Int32), Errors.LineTooLong)
    # Read fields incrementally. `i` is the index to read from, beginning at 1.
    # The parse_x functions will return the next index to begin from (after the tab)
    (qname, i) = @? parse_str(mem, 1)
    (qlen, i) = @? parse_int(mem, i, false)
    (qstart, i) = @? parse_int(mem, i, true)
    # qend may only be zero for unmapped records, e.g. minimap2's `--paf-no-hit`
    # lines, so peek at the strand to check that before parsing it below.
    (qend, i) = @? parse_int(mem, i, true)
    if iszero(qend) && (i > lastindex(mem) || @inbounds(mem[i]) != UInt8('*'))
        return ParserException((i - 1) % Int32, Errors.InvalidZero)
    end
    qend ≥ qstart || return ParserException(i % Int32, Errors.BackwardsIndices)
    qend > qlen && return ParserException(i % Int32, Errors.PositionOutOfBounds)

    # Load the strand field. We need at least 14 more bytes to encode 7 more mandatory fields
    # plus 7 more tabs.
    i + 14 > lastindex(mem) &&
        return ParserException(lastindex(mem) % Int32, Errors.TooFewFields)
    b = @inbounds mem[i]
    strand = if b == UInt8('*')
        return finish_unmapped!(record, mem, qname, qlen, i + 2)
    elseif b == UInt8('-')
        0x02
    elseif b == UInt8('+')
        0x01
    else
        return ParserException(i % Int32, Errors.InvalidStrand)
    end
    @inbounds mem[i + 1] == UInt8('\t') ||
        return ParserException((i + 1) % Int32, Errors.TooFewFields)

    # Load rest of the fields
    (tname, i) = @? parse_str(mem, i + 2) # i + 2 because strand plus tab took two bytes
    (tlen, i) = @? parse_int(mem, i, false)
    (tstart, i) = @? parse_int(mem, i, true)
    (tend, i) = @? parse_int(mem, i, false)
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

    # Fill in fields of the struct. The aux data begins at `i`.
    # Note: PAF uses zero-based semi-open indexing like e.g. Python,
    # so compensate for that in qstart and qend.
    record.data = mem
    record.qname_len = length(qname) % Int32
    record.tname_start = first(tname) % Int32
    record.tname_len = length(tname) % Int32
    record.aux_start = i % Int32
    record.qlen = qlen
    record.qstart = qstart + one(qstart)
    record.qend = qend
    record.tlen = tlen
    record.tstart = tstart + one(tstart)
    record.tend = tend
    record.matches = matches
    record.alnlen = alnlen
    record.mapq = mapq
    record.strand = strand

    return record
end

# If the record is unmapped we don't even bother reading
# the rest of the fields, so take this fast path
# i is the byte index of the strand plus two
function finish_unmapped!(
        record::PAFRecord,
        mem::ImmutableMemoryView{UInt8},
        qname::UnitRange{Int},
        qlen::Int32,
        i::Int,
    )
    # Determine the position of aux data
    for _ in 1:6
        i = findnext(==(UInt8('\t')), mem, i)
        i = if isnothing(i)
            return ParserException(lastindex(mem) % Int32, Errors.TooFewFields)
        else
            i + 1
        end
    end
    i = findnext(==(UInt8('\t')), mem, i)
    aux_start = isnothing(i) ? lastindex(mem) + 1 : i + 1

    # Fill in fields
    record.data = mem
    record.qname_len = length(qname) % Int32
    record.tname_start = 1
    record.tname_len = 0
    record.aux_start = aux_start % Int32
    record.qlen = qlen
    record.qstart = 0
    record.qend = 0
    record.tlen = 0
    record.tstart = 0
    record.tend = 0
    record.matches = 0
    record.alnlen = 0
    record.mapq = 0xff
    record.strand = 0x00

    return record
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
    )::Union{Tuple{Int32, Int}, ParserException}
    # Accumulate in an Int, which can't overflow because we check against
    # typemax(Int32) after every digit.
    n = 0
    i = from
    for outer i in from:lastindex(v)
        b = @inbounds v[i]
        # If we hit a tab, this field is over.
        if b == UInt8('\t')
            i == from && return ParserException(i % Int32, Errors.EmptyInteger)
            (!allow_zero & iszero(n)) &&
                return ParserException(i % Int32, Errors.InvalidZero)
            return (n % Int32, i + 1)
        end
        b ∈ 0x30:0x39 || return ParserException(i % Int32, Errors.BadInteger)
        n = 10 * n + (b - 0x30)
        n > typemax(Int32) && return ParserException(i % Int32, Errors.IntegerOverflow)
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
        (n % Int32, i + 1)
    end
end

"""
    PAFReader(io::Union{IO, AbstractBufReader}; copy::Bool=true)

Construct a `PAFReader`, an iterator of [`PAFRecord`](@ref) read from `io`.
Iterating this object returns a `PAFRecord` for each valid line of input,
and throws a [`ParserException`](@ref) when an invalid record is found.
For efficiency, the auxiliary fields are not validated until they are accessed.

If `io` is an `IO`, it is wrapped in a `BufferIO.BufReader` with a buffer size
given by the keyword `buf_size`.

If `copy` is false, the *same* record will be returned on each iteration.
This record is a view into the reader's buffer, and is only valid until the reader is
advanced. This removes a few allocations and copies per iteration, but iterating the
reader will overwrite old records.

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
           fst = first(reader)
           all(reader) do record
               record === fst
           end
       end
true
```

See also: [`PAFRecord`](@ref), [`try_next!`](@ref)
"""
mutable struct PAFReader{I <: AbstractBufReader}
    const io::ExpandReader{I}
    const record::PAFRecord
    # Number of bytes of the last parsed line, including its newline, which have not
    # yet been consumed from `io`. Consuming is deferred to the next call to `try_next!`,
    # such that the record, which is a view into the buffer of `io`, stays valid.
    to_consume::Int
    line::Int
    # Whether to iterate copies of the record stored in this struct
    copy::Bool
end

function PAFReader(io::AbstractBufReader; copy::Bool = true)
    return PAFReader{typeof(io)}(ExpandReader(io), PAFRecord(), 0, 1, copy)
end

function PAFReader(io::IO; buf_size::Int = 2^16, copy::Bool = true)
    return PAFReader(BufReader(io, max(buf_size, 16)); copy)
end

Base.IteratorSize(::Type{<:PAFReader}) = Base.SizeUnknown()
Base.eltype(::Type{<:PAFReader}) = PAFRecord
Base.close(reader::PAFReader) = close(reader.io) # TODO: Docstring

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
        rec = reader.copy ? copy(res) : res
        (rec, nothing)
    end
end

"""
    try_next!(reader::PAFReader) -> Union{Nothing, PAFRecord, ParserException}

Try to read a record from the [`PAFReader`](@ref) and advance the reader. This may yield one of three options:
* If the reader is empty, return `nothing` and do not advance the reader
* If the next record is invalid, return (do not throw) a `ParserException` and
  do not advance the reader
* If the next record is valid, return it as a [`PAFRecord`](@ref) and advance the reader

Even if the reader itself is not advanced, it may still fill its internal buffer, advancing
the `IO` object it wraps.

The returned record is the reader's own record, which is a view into the reader's buffer.
It is overwritten, and its data invalidated, by the next call to `try_next!`.

# Examples
```jldoctest
julia> reader = PAFReader(open(path_to_paf));

julia> try_next!(reader) isa PAFRecord
true

julia> [try_next!(reader) for i in 1:2] isa Vector{PAFRecord}
true

julia> typeof([try_next!(reader) for i in 1:100])
Vector{Nothing} (alias for Array{Nothing, 1})

julia> close(reader); reader = PAFReader(IOBuffer("bad data"));

julia> try_next!(reader) isa PAF.ParserException
true

julia> all(i -> try_next!(reader) isa PAF.ParserException, 1:100)
true
```
"""
function try_next!(reader::PAFReader)::Union{Nothing, PAFRecord, ParserException}
    io = reader.io
    # Consume the previous line, now that the previous record is invalidated.
    if !iszero(reader.to_consume)
        consume(io, reader.to_consume)
        reader.to_consume = 0
    end
    buffer = get_buffer(io)
    newline = findfirst(==(0x0a), buffer)
    # Keep reading until we find a newline or reach EOF
    while newline === nothing
        # The number of bytes in buffer we have already memchrd searching for a newline
        n_scanned = length(buffer)
        # This signals EOF - we might still have data in the buffer, if the underlying stream
        # ends with a record that doesn't end in a newline
        iszero(fill_buffer(io)) && break
        buffer = get_buffer(io)
        newline = findnext(==(0x0a), buffer, n_scanned + 1)
    end
    has_windows_newline = false
    # Last index before \n, and before \r if the newline is \r\n as on Windows
    last_index = if newline === nothing
        length(buffer)
    else
        if isone(newline)
            0
        else
            has_windows_newline = @inbounds(buffer[newline - 1]) == UInt8('\r')
            newline - (1 + has_windows_newline)
        end
    end
    parse_mem = buffer[1:last_index]
    # If there is no memory to parse (and we have already verified that
    # there are no more bytes in the IO), we end
    # This is only the end of the file if the rest of the stream is whitespace.
    # The buffer may end after the empty line, so keep filling it until EOF.
    if isempty(parse_mem)
        line = reader.line
        col = 0
        n_scanned = 0
        while true
            @inbounds for b in buffer[(n_scanned + 1):end]
                col += 1
                line += b == UInt8('\n')
                col *= b != UInt8('\n')
                if b ∉ (UInt8('\r'), UInt8('\n'), UInt8('\t'), UInt8(' '))
                    return ParserException(Errors.EmptyNonTailingLine, col % Int32, line)
                end
            end
            n_scanned = length(buffer)
            iszero(fill_buffer(io)) && return nothing
            buffer = get_buffer(io)
        end
    end
    # If this results in an error, we return the error.
    res = parse_line!(reader.record, parse_mem)
    res isa ParserException &&
        return ParserException(res.kind, getfield(res, :index_in_line), reader.line)
    # Else, we update the reader. This means `next!` will continuously return
    # an error if an error is found.
    # If the newline was found, we compensate for the newline by also consuming that byte,
    # hence the extra `newline !== nothing`.
    reader.to_consume = last_index + has_windows_newline + (newline !== nothing)
    reader.line += 1
    return res
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
