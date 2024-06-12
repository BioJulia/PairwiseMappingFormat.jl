# TODO: Module doc
module PairwiseMappingFormat

export PAFReader,
    PAFRecord,
    aln_identity,
    aux_data,
    try_next!,
    query_coverage,
    target_coverage

# Not exported, because this very PAF-specific, and operates on foreign types.
public try_parse

using StringViews: StringView
using MemViews: MemView, ImmutableMemView

using SAMAuxData.SAM: Auxiliary

"""
    PAFRecord(buffer_size::Int=2)

Mutable type representing a line PAF line. The content of the record
is accessed through its properties.

# Extended help
The following properties may be used:
* `qname` and `tname` return `StringView`s, and give the name
  if the query and target. These may be empty or contain any
  bytes except a tab
* `qlen` and `tlen` are `Int`s, and gives the length of the query and target
  sequences, respectively. This must be > 0
* `qstart`, `qend`, `tstart`, `tend` are `Int`s, and give the starting and
  ending positions of the alignments on the query and target strand.
  These uses one-based, closed interval indexing like normal in Julia.
  The ending positions are always >= the starting ones
* `matches::Int` gives the number of nucleotides / residues which match
  (i.e. are equal) in the alignment.
* `alnlen::Int` gives the length of the alignment
* `mapq::Union{Int, Nothing}` gives the mapping quality, or `nothing` if 
  this information is unavailable
* `is_rc::Bool` tells if the query and target strands are reverse-complement
  relative to each other
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
    qstart::Int32 # always >= 1
    qend::Int32 # always >= qstart
    tlen::Int32 # always >= 1
    tstart::Int32 # always >= 1
    tend::Int32 # always >= tstart
    matches::Int32
    alnlen::Int32 # always >= 1
    mapq::UInt8 # stored as 0xff for missing
    is_rc::Bool
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
        getfield(record, :is_rc)
    else
        error(lazy"Type PAFRecord has no property $sym")
    end
end

# For tab completion
Base.propertynames(::PAFRecord) = (:qname, :tname, :mapq, :qlen, :qstart, :qend, :tlen, :tstart, :tend, :matches, :alnlen, :is_rc)

function PAFRecord(size::Int=2)
    data = Vector{UInt8}(undef, max(size, 0))
    PAFRecord(data, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0xff, false)
end

function Base.show(io::IO, ::MIME"text/plain", record::PAFRecord)
    buf = IOBuffer()
    println(buf, "PAFRecord:")
    println(buf, "  Query:    ", qname(record))
    println(buf, "  Target:   ", tname(record))
    println(buf, "  Q cov:    ", round(query_coverage(record); digits=4))
    println(buf, "  T cov:    ", round(target_coverage(record); digits=4))
    println(buf, "  Identity: ", round(aln_identity(record); digits=4))
    qual = mapq(record)
    qual = qual === nothing ? "255 (missing)" : string(qual)
    println(buf, "  Quality:  ", qual)
    print(buf, "  Aux data: ")
    write(buf, repr_aux(record))
    print(io, String(take!(buf)))
end

# Print the AUX fields indented
function repr_aux(record::PAFRecord)
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
function tname(record::PAFRecord)
    ql = getfield(record, :qname_len)
    span = (ql+ 1):(ql + getfield(record, :tname_len))
    StringView(ImmutableMemView(getfield(record, :data))[span])
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
"""
function aln_identity(record::PAFRecord)::Float64
    record.matches / record.alnlen
end

"""
    query_coverage(rec::Record) -> Float64

Compute the approximate fraction of the query covered by the alignment.
This is computed as the alignment length divided by the query length,
and thus may be inaccurate if there are deletions in the alignment.
"""
function query_coverage(record::PAFRecord)::FLoat64
    (record.qend - record.qstart + Int32(1)) / record.qlen
end

"""
    target_coverage(rec::Record) -> Float64

Compute the approximate fraction of the target covered by the alignment.
This is computed as the alignment length divided by the target length,
and thus may be inaccurate if there are deletions in the alignment.
"""
function target_coverage(record::PAFRecord)::FLoat64
    (record.tend - record.tstart + Int32(1)) / record.tlen
end

"""
    aux_data(rec::Record) -> SAMAuxData.SAM.Auxiliary
"""
function aux_data(record::PAFRecord)
    Auxiliary(
        getfield(record, :data),
        getfield(record, :qname_len) + getfield(record, :tname_len) + 1,
    )
end

# TODO: Doc
module Errors
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

# TODO: Error doc

# TODO: Doc
struct ParserError
    index_in_line::Int32
    kind::Err
end

# TODO: Doc
struct ParserException <: Exception
    error::ParserError
    line::Int
end

function Base.showerror(io::IO, err::ParserException)
    buf = IOBuffer()
    print(buf, "Error when parsing PAF record on line ")
    print(buf, err.line)
    print(buf, ", near byte number ")
    print(buf, err.error.index_in_line)
    print(buf, " in line: ")
    kind = err.error.kind
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
        res isa ParserError ? (return res) : res
    end
end

# TODO: Doc
try_parse(x) = try_parse(ImmutableMemView(x))
# At least 19 bytes of input are stored in the fixed fields, so
# we need at most 19 fewer bytes in the data field of the vector.
try_parse(x::ImmutableMemView{UInt8}) = parse_line!(PAFRecord(length(x) - 19), x)

function Base.parse(::Type{PAFRecord}, s::Union{String, SubString{String}})
    y = try_parse(s)
    y isa PAFRecord ? y : throw(ParserException(y, 1))
end

function parse_line!(
    record::PAFRecord,
    mem::ImmutableMemView{UInt8},
)::Union{PAFRecord, ParserError}
    if lastindex(mem) > typemax(Int32)
        error("PairwiseMappingFormat.jl can't handle lines longer than 2147483647 bytes")
    end
    # Read fields incrementally. `i` is the index to read from, beginning at 1.
    # The parse_x functions will return the next index to begin from (after the tab)
    (qname, i) = @? parse_str(mem, 1)
    (qlen, i) = @? parse_int(mem, i, false)
    (qstart, i) = @? parse_int(mem, i, true)
    (qend, i) = @? parse_int(mem, i, false)
    qend ≥ qstart || return ParserError(i % Int32, Errors.BackwardsIndices)
    qend > qlen && return ParserError(i % Int32, Errors.PositionOutOfBounds)

    # Load the strand field.
    i > lastindex(mem) - 1 && return ParserError(lastindex(mem) % Int32, Errors.TooFewFields)
    b = @inbounds mem[i]
    b ∈ (UInt8('+'), UInt8('-')) || return ParserError(i % Int32, Errors.InvalidStrand)
    @inbounds mem[i + one(i)] == UInt8('\t') ||
              return ParserError((i + 1) % Int32, Errors.TooFewFields)

    # Load rest of the fields
    (tname, i) = @? parse_str(mem, i + 2) # i + 2 because strand plus tab took two bytes
    (tlen, i) = @? parse_int(mem, i, false)
    (tstart, i) = @? parse_int(mem, i, true)
    (tend, i) = @? parse_int(mem, i, false)
    tend ≥ tstart || return ParserError(i % Int32, Errors.BackwardsIndices)
    tend > tlen && return ParserError(i % Int32, Errors.PositionOutOfBounds)
    (matches, i) = @? parse_int(mem, i, true)
    (alnlen, i) = @? parse_int(mem, i, false)
    (mapq, i) = @? parse_int(mem, i, true, true)
    
    # A missing mapq is encoded as 255 in the PAF format, we store it as such
    mapq = if mapq > 255
        return ParserError(i % Int32, Errors.IntegerOverflow)
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
    unsafe_copyto!(dataview, doff, mem, first(qname), length(qname))
    doff += length(qname)
    unsafe_copyto!(dataview, doff, mem, first(tname), length(tname))
    doff += length(tname)
    iszero(auxlen) || unsafe_copyto!(dataview, doff, mem, i, auxlen)
    

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
    record.is_rc = b == UInt8('-')

    record
end

# Just get the index of the next \t
function parse_str(
    v::ImmutableMemView{UInt8},
    from::Int,
)::Union{Tuple{UnitRange{Int}, Int}, ParserError}
    i = findnext(==(UInt8('\t')), v, from)
    if isnothing(i)
        ParserError(lastindex(v) % Int32, Errors.TooFewFields)
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
)::Union{Tuple{Int32, Int}, ParserError}
    n = Int32(0)
    i = from
    for outer i in from:lastindex(v)
        b = @inbounds v[i]
        # If we hit a tab, this field is over.
        if b == UInt8('\t')
            i == from && return ParserError(i % Int32, Errors.EmptyInteger)
            (!allow_zero & iszero(n)) && return ParserError(i % Int32, Errors.InvalidZero)
            return (n, i + 1)
        end
        b ∈ 0x30:0x39 || return ParserError(i % Int32, Errors.BadInteger)
        n < 0 && return ParserError(i % Int32, Errors.IntegerOverflow)
        n = Int32(10) * n + (b - Int32(48))
    end
    # at_end is if this is the mapq field, which does not need to end with a tab
    if !at_end
        ParserError(lastindex(v) % Int32, Errors.TooFewFields)
    elseif from > lastindex(v)
        ParserError(i % Int32, Errors.EmptyInteger)
    elseif !allow_zero & iszero(n)
        ParserError(i % Int32, Errors.InvalidZero)
    else
        (n, i + 1)
    end
end

# TODO: Doc
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

function PAFReader(io::IO; copy::Bool=true)
    # 2^16 and 512 are reasonable buffer sizes, they are arbitrary
    mem = Memory{UInt8}(undef, 2^16)
    rec = PAFRecord(512)
    PAFReader{typeof(io)}(io, rec, mem, 1, 0, 1, copy)
end

Base.IteratorSize(::Type{<:PAFReader}) = Base.SizeUnknown()
Base.eltype(::Type{<:PAFReader}) = PAFRecord
Base.close(reader::PAFReader) = close(reader.io) # TODO: Docstring

function PAFReader(f::Function, io::IO; copy::Bool=true)
    reader = PAFReader(io; copy)
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
        getfield(record, :is_rc),
    )
end

function Base.iterate(reader::PAFReader, ::Nothing=nothing)
    res = try_next!(reader)
    if res isa ParserError
        throw(ParserException(res, reader.line))
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

# TODO: Doc
function try_next!(reader::PAFReader)::Union{Nothing, PAFRecord, ParserError}
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
    isempty(parse_mem) && return nothing
    # If this results in an error, we return the error.
    parse_result = @? parse_line!(reader.record, parse_mem)
    # Else, we update the reader. This means `next!` will continuously return
    # an error if an error is found.
    # If the newline was found, we compensate for the newline by skipping that byte in
    # the next iteration, hence the extra `newline !== nothing`.
    reader.index = last_index + 1 + has_windows_newline + (newline !== nothing)
    reader.line += 1
    return parse_result
end

end # module PairwiseMappingFormat
