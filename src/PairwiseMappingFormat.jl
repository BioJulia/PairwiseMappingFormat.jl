"""
    module PairwiseMappingFormat

Parse files of the PairwisemAppingFormat (PAF) format.
"""
module PairwiseMappingFormat

export PAFReader,
    PAFRecord,
    aln_identity,
    aux_data,
    try_next!,
    query_coverage,
    target_coverage,
    is_mapped

# Not exported, because this very PAF-specific, and operates on foreign types.
public try_parse, ParserException, Errors, Err

using PrecompileTools: @compile_workload, @setup_workload, @recompile_invalidations

# Note: I don't like to do this. Loading StringViews should cause widespread
# invalidation, but it does.
# This is because Base contain type unstable inference with AbstractString,
# which doesn't cause a lot of issues in Base since there aren't a lot of
# AbstractStrings in Base.
@recompile_invalidations begin
    using StringViews: StringView
end

using MemViews: MemView, ImmutableMemView

using XAMAuxData.SAM: Auxiliary

"""
    PAFRecord(buffer_size::Int)

Mutable type representing a single PAF line. The content of the record
is accessed through its properties.

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
The following properties may be used:
* `qname::StringView`. The query name, May be empty, and contain any bytes
  except `\t` and `\n`.
* `tname::Union{StringView, Nothing}`. Target name. Like `qname`, but is `nothing`
  if and only if the record is unmapped.
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
* `is_rc::Union{Bool,Nothing}` tells if the query and target strands are reverse-complement
  relative to each other. Is `nothing` if the record is unmapped.
"""
mutable struct PAFRecord
    # Data contains query name, subject name, AUX data,
    # without delimiters and with no trailing whitespace.
    const data::Vector{UInt8}
    qname_len::Int32
    tname_len::Int32
    qlen::Int32 # always >= 1
    # Note: qstart/qend (and tstart/tend) are zero-indexed semi-open interval
    # in PAF format, but stored in this struct one-based.
    qstart::Int32 # always >= 1, except when unmapped
    qend::Int32 # always >= qstart
    tlen::Int32 # always >= 1, except when unmapped
    tstart::Int32 # always >= 1, except when unmapped
    tend::Int32 # always >= tstart
    matches::Int32
    alnlen::Int32 # always >= 1, except when unmapped
    mapq::UInt8 # stored as 0xff for missing
    strand::UInt8 # 0x00: *, 0x01: -, 0x02: +
    # 16 bits of padding, which we can use for something else later
end

# Since the interface of the records are defined in terms of its
# properties, we need to be able to abstract over them
function Base.getproperty(record::PAFRecord, sym::Symbol)
    if sym === :qname
        qname(record)
    elseif sym == :tname
        tname(record)
    elseif sym == :mapq
        mapq(record)
    elseif sym in (:qlen, :qstart, :qend, :tlen, :tstart, :tend, :matches, :alnlen)
        Int(getfield(record, sym))
    elseif sym == :is_rc
        is_rc(record)
    else
        error(lazy"Type PAFRecord has no property $sym")
    end
end

"""
    is_mapped(x::PAFRecord) -> Bool

Compute whether the `PAFRecord` is mapped. An unmapped record will have
the properties `is_rc` and `tname` unavailable.
The properties `qname` and `qlen`, and the auxiliary data of an unmapped record
can be relied on, but the remaining properties contain arbitrary data.
Note that some PAF parsers do not handle unmapped records correctly, so be
wary when writing unmapped records.

# Examples
```jldoctes
julia> is_mapped(record)
true
```

See also: [`PAFRecord`](@ref)
"""
is_mapped(x::PAFRecord) = !iszero(getfield(x, :strand))

function is_rc(record::PAFRecord)::Union{Bool, Nothing}
    x = getfield(record, :strand)
    iszero(x) ? nothing : x == 0x02
end

# For tab completion
Base.propertynames(::PAFRecord) = (:qname, :tname, :mapq, :qlen, :qstart, :qend, :tlen, :tstart, :tend, :matches, :alnlen, :is_rc)

function PAFRecord(size::Int=0)
    data = Vector{UInt8}(undef, max(size, 0))
    PAFRecord(data, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0xff, 0x00)
end

function Base.show(io::IO, ::MIME"text/plain", record::PAFRecord)
    buf = IOBuffer()
    is_mapped(record) || print(buf, "Unmapped ")
    println(buf, "PAFRecord:")
    println(buf, "  Query:    \"", qname(record), '"')
    if is_mapped(record)
        println(buf, "  Target:   \"", tname(record), '"')
        println(buf, "  Q cov:    ", round(query_coverage(record); digits=4))
        println(buf, "  T cov:    ", round(target_coverage(record); digits=4))
        println(buf, "  Identity: ", round(aln_identity(record); digits=4))
        qual = mapq(record)
        qual = qual === nothing ? "255 (missing)" : string(qual)
        println(buf, "  Quality:  ", qual)
    end
    print(buf, "  Aux data: ")
    write(buf, repr_aux(record))
    print(io, String(take!(buf)))
end

# Print the AUX fields indented
function repr_aux(record::PAFRecord)::Vector{UInt8}
    buf = IOBuffer()
    out = IOBuffer()
    show(buf, MIME"text/plain"(), aux_data(record))
    seekstart(buf)
    for line in eachline(buf; keep=true)
        write(out, "  ", line)
    end
    take!(out)[3:end]
end

# Return StringView to not allocate
function qname(record::PAFRecord)
    StringView(ImmutableMemView(getfield(record, :data))[1:(getfield(record, :qname_len))])
end

# Return StringView to not allocate
function tname(record::PAFRecord)::Union{StringView, Nothing}
    if is_mapped(record)
        ql = getfield(record, :qname_len)
        span = (ql+ 1):(ql + getfield(record, :tname_len))
        StringView(ImmutableMemView(getfield(record, :data))[span])
    else
        nothing
    end
end

# Note: We return nothing instead of missing because fuck missing
function mapq(record::PAFRecord)::Union{Int, Nothing}
    x = getfield(record, :mapq)
    x == 0xff ? nothing : Int(x)
end

"""
    aln_identity(rec::Record) -> Float64

Return the nucleotide / residue identity, defined as the number of
matches divided by the alignment length.
This is a number between 0 and 1.

# Examples
```jldoctest
julia> aln_identity(record)
1.0
```

See also: [`query_coverage`](@ref), [`target_coverage`](@ref)
"""
function aln_identity(record::PAFRecord)::Float64
    record.matches / record.alnlen
end

"""
    query_coverage(rec::Record) -> Float64

Compute the fraction of the query covered by the alignment.

# Examples
```jldoctest
julia> query_coverage(record)
0.9999535124653004
```

See also: [`target_coverage`](@ref), [`aln_identity`](@ref)
"""
function query_coverage(record::PAFRecord)::Float64
    (record.qend - record.qstart + Int32(1)) / record.qlen
end

"""
    target_coverage(rec::Record) -> Float64

Compute the fraction of the target covered by the alignment.

# Examples
```jldoctest
julia> target_coverage(record)
0.044934629307437725
```

See also: [`query_coverage`](@ref), [`aln_identity`](@ref)
"""
function target_coverage(record::PAFRecord)::Float64
    (record.tend - record.tstart + Int32(1)) / record.tlen
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

julia> aux["k1"] = "add a new aux string like this";

julia> haskey(aux, "k1")
true
```
"""
function aux_data(record::PAFRecord)
    Auxiliary(
        getfield(record, :data),
        getfield(record, :qname_len) + getfield(record, :tname_len) + 1,
    )
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

See also: [`ParserException`](@ref)
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
    if s ∈ (:kind, :line)
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
    end
    print(buf, s)
    print(io, String(take!(buf)))
end

# Use this macro to return early and reduce code duplication.
# It's similar to `a = @something b return nothing`
macro var"?"(expr)
    quote
        local res = $(esc(expr))
        res isa ParserException ? (return res) : res
    end
end

"""
    try_parse(x) -> Union{PAFRecord, ParserException}

Try parsing `x` into a `PAFRecord`. `x` may be any type that implements
`MemView`, such as `String` or `Vector{UInt8}`.

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
try_parse(x) = try_parse(ImmutableMemView(x))
# At least 19 bytes of input are stored in the fixed fields, so
# we need at most 19 fewer bytes in the data field of the vector.
try_parse(x::ImmutableMemView{UInt8}) = parse_line!(PAFRecord(length(x) - 19), x)

function Base.parse(::Type{PAFRecord}, s::Union{String, SubString{String}})
    y = try_parse(s)
    y isa PAFRecord ? y : throw(y)
end

function parse_line!(
    record::PAFRecord,
    mem::ImmutableMemView{UInt8},
)::Union{PAFRecord, ParserException}
    if lastindex(mem) > typemax(Int32)
        error("PairwiseMappingFormat.jl can't handle lines longer than 2147483647 bytes")
    end
    # Read fields incrementally. `i` is the index to read from, beginning at 1.
    # The parse_x functions will return the next index to begin from (after the tab)
    (qname, i) = @? parse_str(mem, 1)
    (qlen, i) = @? parse_int(mem, i, false)
    (qstart, i) = @? parse_int(mem, i, true)
    (qend, i) = @? parse_int(mem, i, false)
    qend ≥ qstart || return ParserException(i % Int32, Errors.BackwardsIndices)
    qend > qlen && return ParserException(i % Int32, Errors.PositionOutOfBounds)

    # Load the strand field. We need at least 14 more bytes to encode 7 more mandatory fields
    # plus 7 more tabs.
    i + 14 > lastindex(mem) && return ParserException(lastindex(mem) % Int32, Errors.TooFewFields)
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
    (alnlen, i) = @? parse_int(mem, i, false)
    (mapq, i) = @? parse_int(mem, i, true, true)
    
    # A missing mapq is encoded as 255 in the PAF format, we store it as such
    mapq = if mapq > 255
        return ParserException(i % Int32, Errors.IntegerOverflow)
    else
        mapq % UInt8
    end

    # Copy data to the .data field - query name, target name, AUX data
    auxlen = lastindex(mem) - i + 1
    filled = length(qname) + length(tname) + auxlen
    data = getfield(record, :data)
    length(data) == filled || resize!(data, filled)
    doff = 1
    dataview = MemView(data)
    unsafe_copyto!(dataview, mem[qname])
    doff += length(qname)
    unsafe_copyto!(dataview[doff:end], mem[tname])
    doff += length(tname)
    iszero(auxlen) || unsafe_copyto!(dataview[doff:end], mem[i:end])
    

    # Fill in fields of the struct
    # Note: PAF uses zero-based semi-open indexing like e.g. Python,
    # so compensate for that in qstart and qend.
    record.qname_len = length(qname) % Int32
    record.tname_len = length(tname) % Int32
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

    record
end

# If the record is unmapped we don't even bother reading
# the rest of the fields, so take this fast path
# i is the byte index of the strand plus two
function finish_unmapped!(
    record::PAFRecord,
    mem::ImmutableMemView{UInt8},
    qname::UnitRange{Int},
    qlen::Int32,
    i::Int
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
    aux_end = lastindex(mem)
    aux_start = isnothing(i) ? aux_end + 1 : i + 1
    aux = aux_start:aux_end

    # Copy data
    filled = length(qname) + length(aux)
    data = getfield(record, :data)
    dataview = MemView(data)
    length(data) == filled || resize!(data, filled)
    unsafe_copyto!(dataview, mem[qname])
    unsafe_copyto!(dataview[length(qname)+1:end], mem[aux])

    # Fill in fields
    record.qname_len = length(qname) % Int32
    record.tname_len = 0
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

    record

end

# Just get the index of the next \t
function parse_str(
    v::ImmutableMemView{UInt8},
    from::Int,
)::Union{Tuple{UnitRange{Int}, Int}, ParserException}
    i = findnext(==(UInt8('\t')), v, from)
    if isnothing(i)
        ParserException(lastindex(v) % Int32, Errors.TooFewFields)
    else
        (from:i-1, i + 1)
    end
end

# TODO: 40% of runtime is spent in this function. Microoptimize it
function parse_int(
    v::ImmutableMemView{UInt8},
    from::Int,
    allow_zero::Bool,
    at_end::Bool=false,
)::Union{Tuple{Int32, Int}, ParserException}
    n = Int32(0)
    i = from
    for outer i in from:lastindex(v)
        b = @inbounds v[i]
        # If we hit a tab, this field is over.
        if b == UInt8('\t')
            i == from && return ParserException(i % Int32, Errors.EmptyInteger)
            (!allow_zero & iszero(n)) && return ParserException(i % Int32, Errors.InvalidZero)
            return (n, i + 1)
        end
        b ∈ 0x30:0x39 || return ParserException(i % Int32, Errors.BadInteger)
        n < 0 && return ParserException(i % Int32, Errors.IntegerOverflow)
        n = Int32(10) * n + (b - Int32(48))
    end
    # Note: It is not possible to get an InvalidZero error here, because
    # the last field (mapq) does allow zero, and for any other fields,
    # if they do not end with a \t and thus reach this line, a TooFewFields
    # error will be returned.
    # at_end is if this is the mapq field, which does not need to end with a tab
    if !at_end
        ParserException(lastindex(v) % Int32, Errors.TooFewFields)
    elseif from > lastindex(v)
        ParserException(i % Int32, Errors.EmptyInteger)
    else
        (n, i + 1)
    end
end

"""
    PAFReader(io::IO; copy::Bool=true)

Construct a `PAFReader`, an iterator of [`PAFRecord`](@ref) read from `io`.
Iterating this object returns a `PAFRecord` for each valid line of input,
and throws a [`ParserException`](@ref) when an invalid record is found.
For efficiency, the auxiliary fields are not validated until they are accessed.

If `copy` is false, the *same* record will be returned on each iteration,
with its content being overwritten. This removes a few allocations per iteration,
but may cause problems if records of old iterations are stored.

# Examples
```jldoctest
julia> reader = PAFReader(open(path_to_paf)); typeof(reader)
PAFReader{IOStream}

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
mutable struct PAFReader{I <: IO}
    const io::I
    const record::PAFRecord
    mem::Memory{UInt8}
    # Index of first unused byte in buffer. Equal to filled + 1
    # if the buffer has no usable data.
    index::Int32
    # mem[1:filled] contain usable data
    filled::Int32
    line::Int
    # Whether to iterate copies of the record stored in this struct
    copy::Bool
end

function PAFReader(io::IO; buf_size::Int=2^16, copy::Bool=true)
    mem = Memory{UInt8}(undef, max(buf_size, 16))
    rec = PAFRecord(512)
    PAFReader{typeof(io)}(io, rec, mem, 1, 0, 1, copy)
end

Base.IteratorSize(::Type{<:PAFReader}) = Base.SizeUnknown()
Base.eltype(::Type{<:PAFReader}) = PAFRecord
Base.close(reader::PAFReader) = close(reader.io) # TODO: Docstring

function PAFReader(f::Function, io::IO; kwargs...)
    reader = PAFReader(io; kwargs...)
    try
        f(reader)
    finally
        close(reader)
    end
end

function Base.copy(record::PAFRecord)
    PAFRecord(
        copy(getfield(record, :data)),
        getfield(record, :qname_len),
        getfield(record, :tname_len),
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

function Base.iterate(reader::PAFReader, ::Nothing=nothing)
    res = try_next!(reader)
    if res isa ParserException
        throw(res)
    elseif res === nothing
        nothing
    else
        rec = reader.copy ? copy(res) : res
        (rec, nothing)
    end
end

# This is never called during normal operations, and only if a single line
# does not fit in the default buffer, which is highly unlikely.
function reallocate!(reader::PAFReader)
    mem = reader.mem
    newlen = length(mem) * 2
    if newlen >= typemax(Int32)
        error("PAFReader encountered a line with more than 2147483646 bytes.")
    end
    new = Memory{UInt8}(undef, newlen)
    unsafe_copyto!(new, 1, mem, 1, reader.filled)
    reader.mem = new
    reader
end

# Shift the data inside the buffer such that the usable
# data begins at index 1, thus making room for more
# data at the end of the buffer
function shift_buffer!(reader::PAFReader)
    isone(reader.index) && return reader
    mem = reader.mem
    len = reader.filled - reader.index + 1
    copyto!(mem, 1, mem, reader.index, len)
    reader.index = 1
    reader.filled = len
    reader
end

# Read data into the buffer. Returns zero if and only if underlying
# stream is EOF (or buggy).
function read_buffer!(reader::PAFReader)::Int
    mem = reader.mem
    remaining = length(mem) - reader.filled
    # If no space, try to shift data in the buffer to make room
    if iszero(remaining)
        shift_buffer!(reader)
        remaining = length(mem) - reader.filled
        # If still no space, the buffer is truly full, and we need to 
        # resize it
        if iszero(remaining)
            reallocate!(reader)
            mem = reader.mem # it was change in reallocate!
            remaining = length(mem) - reader.filled
        end
    end
    # This part annoys me a lot. Julia has no efficient API to read bytes
    # into a byte buffer. How on Earth can that be true?
    # * read! will throw EOF if the memory is too short, and there is no
    #   way to query how many bytes are available
    # * readbytes! is not optimised for Memory
    # * unsafe_read has almost no promised semantics so practically can't be used.
    # The approach here uses unsafe_read in a manner similar to TranscodingStreams.jl.
    # It uses bytesavailable to check how many bytes to read and avoid buffer overflows etc,
    # and when there are zero bytes available, it reads a single byte, hoping that
    # more bytes will be available.
    if !iszero(remaining)
        n_available = bytesavailable(reader.io)
        # Read one byte and hope more bytes become available
        if iszero(n_available) && !eof(reader.io)
            byte = read(reader.io, UInt8)
            @inbounds mem[reader.filled + 1] = byte
            reader.filled += 1
            remaining -= 1
            n_available = bytesavailable(reader.io)
        end
        # Read with unsafe_read
        n_bytes = min(n_available, remaining)
        if !iszero(n_bytes)
            GC.@preserve mem unsafe_read(reader.io, pointer(mem, reader.filled + 1), n_bytes)
            reader.filled += n_bytes
        end
        return n_bytes
    else
        return 0
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
    # n_scanned: The number of bytes in buffer we have already memchrd searching for a newline
    n_scanned = 0
    newline = findfirst(
        ==(0x0a),
        ImmutableMemView(reader.mem)[(reader.index + n_scanned):(reader.filled % Int)],
    )
    n_bytes_read = 0
    # Keep reading until we find a newline or reach EOF
    while newline === nothing
        n_scanned += reader.filled - reader.index - n_scanned + 1
        n_bytes_read = read_buffer!(reader)
        # This signals EOF - we might still have data in the buffer, if the underlying stream
        # ends with a record that doesn't end in a newline
        iszero(n_bytes_read) && break
        newline = findfirst(
            ==(0x0a),
            ImmutableMemView(reader.mem)[(reader.index + n_scanned):(reader.filled % Int)],
        )
    end
    # If a newline was found, the current value of newline was relative to where
    # the last search began. So, update it to account for the bytes searched in
    # previous iterations of search for tab
    newline === nothing || (newline += reader.index + n_scanned - 1)
    mem = reader.mem
    has_windows_newline = false
    # Last index before \n, and before \r if the newline is \r\n as on Windows
    last_index = if newline === nothing
        reader.filled
    else
        if isone(newline)
            0
        else
            has_windows_newline = @inbounds(mem[newline - 1]) == UInt8('\r')
            newline - (1 + has_windows_newline)
        end
    end
    parse_mem = ImmutableMemView(mem)[reader.index:last_index]
    # If there is no memory to parse (and we have already verified that
    # there are no more bytes in the IO), we end
    if isempty(parse_mem)
        line = reader.line
        col = 0
        @inbounds for b in ImmutableMemView(mem)[reader.index:reader.filled]
            col += 1
            line += b == UInt8('\n')
            col *= b != UInt8('\n')
            if b ∉ (UInt8('\r'), UInt8('\n'), UInt8('\t'), UInt8(' '))
                return ParserException(Errors.EmptyNonTailingLine, col % Int32, line)
            end
        end
        return nothing
    end
    # If this results in an error, we return the error.
    res = parse_line!(reader.record, parse_mem)
    res isa ParserException && return ParserException(res.kind, getfield(res, :index_in_line), reader.line)
    # Else, we update the reader. This means `next!` will continuously return
    # an error if an error is found.
    # If the newline was found, we compensate for the newline by skipping that byte in
    # the next iteration, hence the extra `newline !== nothing`.
    reader.index = last_index + 1 + has_windows_newline + (newline !== nothing)
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
