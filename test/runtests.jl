module TestPAF

using PairwiseMappingFormat
using PairwiseMappingFormat: Errors, try_parse
using Test
using FormatSpecimens
using MemoryViews
using Aqua
using BufferIO: CursorReader, AbstractBufReader
import BufferIO: get_buffer, fill_buffer, consume

# A reader that delivers data in fixed chunks, like a pipe would
mutable struct ChunkReader <: AbstractBufReader
    chunks::Vector{Vector{UInt8}}
    buf::Vector{UInt8}
end
ChunkReader(chunks) = ChunkReader(map(c -> collect(codeunits(c)), chunks), UInt8[])
get_buffer(r::ChunkReader) = ImmutableMemoryView(r.buf)
function fill_buffer(r::ChunkReader)
    isempty(r.chunks) && return 0
    c = popfirst!(r.chunks)
    append!(r.buf, c)
    return length(c)
end
consume(r::ChunkReader, n::Int) = (deleteat!(r.buf, 1:n); nothing)
Base.close(::ChunkReader) = nothing

using XAMAuxData: XAMAuxData
using XAMAuxData.SAM: AuxTag

function cmp_aux(record::PAFRecord, other)
    aux = sort!(collect(aux_data(record)); by = first)
    other = sort!(collect(other); by = first)
    @test length(aux) == length(other)
    for ((kaux, vaux), (ko, vo)) in zip(aux, other)
        @test kaux == AuxTag(ko)
        if vaux isa Union{Real, AbstractString, Char, XAMAuxData.Error}
            @test vaux == vo
        elseif vaux isa AbstractVector
            @test eltype(vaux) == eltype(vo)
            @test vaux == vo
        else
            @warn typeof(vaux)
            @test false
        end
    end
    return
end

@testset "Normal records" begin
    (rec1, rec2, rec3) = PAFReader(collect, open(joinpath(path_of_format("PAF"), "good1.paf")))

    @testset "Field access" begin
        @test rec1.qname == "NODE_1_length_301156_cov_79.148382"
        @test rec1.qlen == 301156
        @test rec1.qstart == 4 # one-based indexing
        @test rec1.is_rc == false
        @test rec1.tname == "CP004047.1"
        @test rec1.tlen == 6701780
        @test rec1.tstart == 2764861 # one-based indexing
        @test rec1.tend == 3066002
        @test rec1.matches == 301142
        @test rec1.alnlen == 301142
        @test rec1.mapq == 0

        # tname and is_rc are type stable, since they throw for unmapped records
        @test (@inferred (r -> r.tname)(rec1)) isa AbstractString
        @test (@inferred (r -> r.is_rc)(rec1)) isa Bool

        @test_throws Exception rec1.foobar
        @test_throws Exception rec1.strand
        @test_throws Exception rec1.data
        @test_throws Exception rec1.qname_len
        @test_throws Exception rec1.tname_len
    end

    @testset "Aux data" begin
        cmp_aux(rec1, ["tp" => 'P', "cm" => 29990, "s1" => 301142, "s2" => 301142, "dv" => 0, "rl" => 0])
        cmp_aux(rec2, ["tp" => 'P', "cm" => 29957, "s1" => 299269, "s2" => 299269, "dv" => 0, "rl" => 38])
        cmp_aux(rec3, ["tp" => 'P', "cm" => 28871, "s1" => 288618, "s2" => 288618, "dv" => 0, "rl" => 57])
    end
end

@testset "Coverage and identity" begin
    function get_rec(qlen::Int, tlen::Int, matches::Int, alnlen::Int)
        fields = split(
            "my_qname\t301156\t3\t11\t+\tmy_tname\t6701780\t7\t13\t301142\t301142\t0",
            '\t'
        )
        fields[2] = string(qlen)
        fields[7] = string(tlen)
        fields[10] = string(matches)
        fields[11] = string(alnlen)
        parse(PAFRecord, join(fields, '\t'))
    end
    @test isapprox(aln_identity(get_rec(500, 400, 325, 340)), 325 / 340)
    @test isapprox(aln_identity(get_rec(12, 16, 9, 11)), 9 / 11)
    @test isapprox(aln_identity(get_rec(100, 1000, 0, 99)), 0)

    @test isapprox(query_coverage(get_rec(511, 398, 325, 340)), 8 / 511)
    @test isapprox(query_coverage(get_rec(61, 19, 9, 9)), 8 / 61)
    @test isapprox(query_coverage(get_rec(17, 16, 15, 16)), 8 / 17)

    @test isapprox(target_coverage(get_rec(511, 398, 325, 340)), 6 / 398)
    @test isapprox(target_coverage(get_rec(61, 19, 9, 9)), 6 / 19)
    @test isapprox(target_coverage(get_rec(17, 16, 15, 16)), 6 / 16)
end

function with_replaced_str(s::String, field::Integer, new::String)
    fields = split(s, '\t')
    fields[field] = new
    return join(fields, "\t")
end

function with_replaced(s::String, field::Integer, new::String)
    y = try_parse(with_replaced_str(s, field, new))
    return y isa PairwiseMappingFormat.ParserException ? y.kind : y
end

@testset "Edge case records" begin
    (rec1, rec2, rec3) = PAFReader(open(joinpath(path_of_format("PAF"), "good2.paf"))) do reader
        y = first(reader, 3)
        @test isempty(reader)
        y
    end
    @test rec1.qname == "中国科学院大学"
    @test rec1.qlen == typemax(Int32)
    # XAMAuxData 0.2 rejects lowercase hex digits
    xX = pkgversion(XAMAuxData) < v"0.2" ? [0xea, 0x19, 0xfe, 0x1b] : XAMAuxData.Errors.InvalidHex
    cmp_aux(rec1, ["tp" => 'S', "cm" => 29990, "s1" => 301142, "dv" => -13.2442f-18, "rl" => 0, "xX" => xX, "K1" => Int8[15, 9, 22, 127, -12, -128, 0]])
    cmp_aux(rec2, ["dv" => 1.243f-4, "rl" => 38, "kk" => "!=[]! And then~~  "])
    cmp_aux(rec3, ["dv" => +17.0f21, "rl" => 38])

    # Empty record
    rec = PAFReader(first, open(joinpath(path_of_format("PAF"), "good3.paf")))
    @test !is_mapped(rec)
    for sym in (:tname, :is_rc, :mapq, :qstart, :qend, :tlen, :tstart, :tend, :matches, :alnlen)
        @test_throws PairwiseMappingFormat.UnmappedError(sym) getproperty(rec, sym)
    end
    for f in (query_coverage, target_coverage, aln_identity)
        @test_throws PairwiseMappingFormat.UnmappedError(Symbol(f)) f(rec)
    end
    @test rec.qname == "query"
    @test rec.qlen == 56
    @test Dict(aux_data(rec)) == Dict{AuxTag, Any}(
        "tp" => 'S',
        "cm" => 4,
        "s1" => 49,
        "dv" => Float32(0.0117),
        "rl" => 0,
    )

    # Unmapped record as written by minimap2 --paf-no-hit, with zeros
    rec = parse(PAFRecord, "read1\t1523\t0\t0\t*\t*\t0\t0\t0\t0\t0\t0")
    @test !is_mapped(rec)
    @test rec.qname == "read1"
    @test rec.qlen == 1523

    # Unmapped record with too few fieds
    @test try_parse("query\t123\t1\t100\t*\t*\t5\t19").kind == Errors.TooFewFields

    # The fields of unmapped records are validated like for mapped records,
    # except that integers may be zero
    unmapped = "read1\t1523\t0\t0\t*\t*\t0\t0\t0\t0\t0\t0"
    for field in 6:12
        @test with_replaced(unmapped, field, "0") isa PAFRecord
    end
    # Use a longer target name, so emptying a field does not make the line too short
    long_unmapped = with_replaced_str(unmapped, 6, "some_target_name")
    for field in 7:12
        @test with_replaced(long_unmapped, field, "") == Errors.EmptyInteger
        @test with_replaced(unmapped, field, "-1") == Errors.BadInteger
        @test with_replaced(unmapped, field, "abc") == Errors.BadInteger
        @test with_replaced(unmapped, field, "4294967296") == Errors.IntegerOverflow
    end
    @test with_replaced(unmapped, 12, "255") isa PAFRecord
    @test with_replaced(unmapped, 12, "256") == Errors.IntegerOverflow
    @test with_replaced(unmapped, 9, "1") == Errors.PositionOutOfBounds
    @test try_parse("read1\t1523\t0\t0\t*\t*\t10\t5\t4\t0\t0\t0").kind == Errors.BackwardsIndices
    @test try_parse("read1\t1523\t0\t0\t*\t*\t10\t0\t5\t0\t0\t0") isa PAFRecord
    @test try_parse("read1\t1523\t0\t0\t*x*\t10\t0\t5\t0\t0\t0").kind == Errors.TooFewFields
    # Aux data is found after the mapq field
    rec = try_parse("read1\t1523\t0\t0\t*\t*\t0\t0\t0\t0\t0\t0\tcm:i:4")
    @test Dict(aux_data(rec)) == Dict{AuxTag, Any}("cm" => 4)
end

@testset "Errors" begin
    good = "my_qname\t301156\t3\t301145\t+\tmy_tname\t6701780\t2764860\t3066002\t301142\t301142\t0"
    @test with_replaced(good, 1, "my other name") isa PAFRecord

    # Negative numbers, zero numbers
    @test with_replaced(good, 2, "-1") == Errors.BadInteger

    for allow_zero in [3, 8, 10, 11, 12]
        @test with_replaced(good, allow_zero, "0") isa PAFRecord
    end

    for not_allow_zero in [2, 4, 7, 9]
        @test with_replaced(good, not_allow_zero, "0") == Errors.InvalidZero
    end

    # qend and tend must be >= qstart tstart
    @test with_replaced(good, 3, "1000000") == Errors.BackwardsIndices
    @test with_replaced(good, 4, "2") == Errors.BackwardsIndices
    @test with_replaced(good, 8, "3500000") == Errors.BackwardsIndices
    @test with_replaced(good, 9, "2000000") == Errors.BackwardsIndices
    @test with_replaced(good, 3, "200000") isa PAFRecord

    # qend and tend must not > qlen tlen
    @test with_replaced(good, 2, "300000") == Errors.PositionOutOfBounds
    @test with_replaced(good, 4, "310000") == Errors.PositionOutOfBounds
    @test with_replaced(good, 7, "3060002") == Errors.PositionOutOfBounds
    @test with_replaced(good, 9, "6711780") == Errors.PositionOutOfBounds

    # Char not valid
    @test with_replaced(good, 5, "/") == Errors.InvalidStrand

    # Too few columns
    for i in 0:11
        s = join(split(good, '\t')[1:i], '\t')
        y = try_parse(s)
        @test y isa PairwiseMappingFormat.ParserException && y.kind == Errors.TooFewFields
    end

    # Missing field
    for int_field in [2, 3, 4, 7, 8, 9, 10, 11, 12]
        @test with_replaced(good, int_field, "") == Errors.EmptyInteger
    end

    # Mapq too high
    @test with_replaced(good, 12, "256") == Errors.IntegerOverflow

    # Integer overflow, including on the last digit
    @test with_replaced(good, 10, "4294967295").matches == typemax(UInt32)
    @test with_replaced(good, 10, "4294967296") == Errors.IntegerOverflow
    @test with_replaced(good, 10, "9999999999") == Errors.IntegerOverflow
    @test with_replaced(good, 10, "99999999999999999999") == Errors.IntegerOverflow
    @test with_replaced(good, 12, "4294967296") == Errors.IntegerOverflow
    @test with_replaced(good, 2, "4294967296") == Errors.IntegerOverflow

    # Sequences longer than typemax(Int32), up to typemax(UInt32)
    m = typemax(UInt32)
    rec = try_parse("q\t$m\t$(m - 1)\t$m\t+\tt\t$m\t0\t$m\t1\t1\t0")
    @test (rec.qlen, rec.qstart, rec.qend) == (m, m, m)
    @test (rec.tlen, rec.tstart, rec.tend) == (m, 1, m)
    @test query_coverage(rec) == 1 / m
    @test target_coverage(rec) == 1.0

    # Empty interval at the very end of the sequence doesn't overflow
    rec = try_parse("q\t$m\t$m\t$m\t+\tt\t5\t0\t5\t1\t1\t0")
    @test (rec.qstart, rec.qend) == (Int(m) + 1, m)
end

@testset "Bad reader" begin
    # Iterating a reader with bad data will throw
    good = "my_qname\t301156\t3\t301145\t+\tmy_tname\t6701780\t2764860\t3066002\t301142\t301142\t0"

    @test only(PAFReader(collect, IOBuffer(good))) isa PAFRecord
    bad = good[1:(end - 2)]
    @test_throws PairwiseMappingFormat.ParserException only(PAFReader(collect, IOBuffer(bad)))
end

@testset "Reader buffer size" begin
    # Large data where buffer needs to be shifted
    good = "my_qname\t301156\t3\t301145\t+\tmy_tname\t6701780\t2764860\t3066002\t301142\t301142\t0"
    large = join([good for _ in 1:1000], '\n')
    PAFReader(IOBuffer(large); buf_size = 16) do reader
        (allgood, n) = foldl(reader; init = (true, 0)) do (isgood, n), record
            isgood &= record.qlen == 301156
            isgood &= record.mapq == 0
            isgood &= record.tstart == 2764861
            (isgood, n + 1)
        end
        @test n == 1000
        @test allgood
    end
end

@testset "AbstractBufReader" begin
    good = "my_qname\t301156\t3\t301145\t+\tmy_tname\t6701780\t2764860\t3066002\t301142\t301142\t0"
    other = replace(good, "my_qname" => "other")
    data = join([good, other, good], '\n')

    # Records iterated with copy=true do not alias the reader's buffer
    records = PAFReader(collect, CursorReader(data))
    @test [r.qname for r in records] == ["my_qname", "other", "my_qname"]

    # With copy=false, each record is a distinct view into the reader's buffer.
    # CursorReader never moves its data, so the views stay valid here.
    PAFReader(CursorReader(data); copy = false) do reader
        rec = try_next!(reader)
        @test rec.qname == "my_qname"
        rec2 = try_next!(reader)
        @test rec2 !== rec
        @test rec2.qname == "other"
    end
    records = PAFReader(collect, CursorReader(data); copy = false)
    @test [r.qname for r in records] == ["my_qname", "other", "my_qname"]

    # try_next! respects copy: With copy=true, records stay valid after advancing,
    # even when the reader's buffer is overwritten
    reader = PAFReader(ChunkReader([good * "\n", other * "\n"]))
    rec = try_next!(reader)
    @test try_next!(reader).qname == "other"
    @test rec.qname == "my_qname"

    # Parsing a vector aliases it, but copying the record does not
    v = collect(codeunits(good))
    rec = try_parse(v)
    copied = copy(rec)
    v[1:2] .= UInt8('X')
    @test rec.qname == "XX_qname"
    @test copied.qname == "my_qname"
end

# We assume that aux data is correctly tested in its own package, here we just test
# that the aux data has the right span and byte content.
function test_aux_content(record_data::String, target_aux::String)
    rec = parse(PAFRecord, record_data)
    return @test collect(MemoryView(aux_data(rec))) == codeunits(target_aux)
end

@testset "Auxiliary data" begin
    good = "my_qname	301156	3	301145	+	my_tname	6701780	2764860	3066002	301142	301142	0	tp:A:P	kL:Z:325 1	s1:i:301142	x3:H:a40140	dv:f:0.0000	rl:i:0"
    goodaux = "tp:A:P	kL:Z:325 1	s1:i:301142	x3:H:a40140	dv:f:0.0000	rl:i:0"
    test_aux_content(good, goodaux)

    stringaux_with_whitespace = "\tHA:Z:lsjd   "
    test_aux_content(good * stringaux_with_whitespace, goodaux * stringaux_with_whitespace)
end

@testset "Parsing lines" begin
    good = "my_qname\t301156\t3\t301145\t+\tmy_tname\t6701780\t2764860\t3066002\t301142\t301142\t0"
    expected = parse(PAFRecord, good)
    same(rec) = rec isa PAFRecord && getfield(rec, :data) == getfield(expected, :data)

    # parse accepts the same inputs as try_parse
    for x in (good, SubString(good, 1), codeunits(good), collect(codeunits(good)), MemoryView(good))
        @test same(parse(PAFRecord, x))
        @test same(try_parse(x))
    end
    @test_throws PairwiseMappingFormat.ParserException parse(PAFRecord, codeunits("abc"))

    # One optional trailing newline is stripped, also before aux data
    @test same(try_parse(good * "\n"))
    @test same(try_parse(good * "\r\n"))
    @test aux_data(parse(PAFRecord, good * "\tab:Z:x\r\n"))["ab"] == "x"

    # Any other newline is an error
    @test try_parse(good * "\n\n").kind == Errors.EmbeddedNewline
    @test try_parse(good * "\r\n\n").kind == Errors.EmbeddedNewline
    @test try_parse(good * "\n" * good).kind == Errors.EmbeddedNewline
    @test try_parse(replace(good, "my_qname" => "my\nqname")).kind == Errors.EmbeddedNewline
    @test try_parse(good * "\tab:Z:x\ny").kind == Errors.EmbeddedNewline
    @test try_parse("\n").kind == Errors.TooFewFields

    # An embedded \r is data
    @test parse(PAFRecord, replace(good, "my_qname" => "my\rqname")).qname == "my\rqname"

    # Other trailing whitespace is not stripped, but part of the last field
    @test try_parse(good * " ").kind == Errors.BadInteger
    @test try_parse(good * "\r").kind == Errors.BadInteger
    # A trailing tab begins an empty aux field, which is only invalid when accessed
    @test try_parse(good * "\t") isa PAFRecord
    @test try_parse(good * " \n").kind == Errors.BadInteger
    @test aux_data(parse(PAFRecord, good * "\tab:Z:x \r\n"))["ab"] == "x "
    # A trailing \r is kept in aux data, where it is invalid for Z values
    @test collect(MemoryView(aux_data(parse(PAFRecord, good * "\tab:Z:x\r")))) == codeunits("ab:Z:x\r")
end

# Note: We don't trim whitespace off ends of records, because whitespace is allowed in
# string AuxData. Trailing whitespace is part of the last field, and may make it invalid.
# After the records, a file may end with whitespace-only lines.
@testset "Trailing whitespace" begin
    data = "my_qname	301156	3	301145	+	my_tname	6701780	2764860	3066002	301142	301142	0	an:Z: dljdl \r\n"
    data *= "some	5	3	4	+	other	9	5	8	2	2	0\r\n"

    reader = PAFReader(IOBuffer(data))
    (rec1, rec2) = collect(reader)

    @test aux_data(rec1)["an"] == " dljdl "

    data *= "\n\n"
    @test !isempty(PAFReader(collect, IOBuffer(data)))

    data *= "\nV\n"
    PAFReader(IOBuffer(data)) do reader
        @test first(reader, 2) isa Vector{PAFRecord}
        @test try_next!(reader).kind == Errors.EmptyNonTailingLine
    end

    # Blank line in the middle, where the buffer ends right after the blank line
    good = "my_qname\t301156\t3\t301145\t+\tmy_tname\t6701780\t2764860\t3066002\t301142\t301142\t0"
    reader = PAFReader(ChunkReader([good * "\n\n", good * "\n" * good * "\n"]))
    @test try_next!(reader) isa PAFRecord
    @test try_next!(reader).kind == Errors.EmptyNonTailingLine

    # Trailing whitespace delivered in multiple chunks is still fine
    reader = PAFReader(ChunkReader([good * "\n\n", "\r\n", "\n"]))
    @test length(collect(reader)) == 1

    # The first whitespace-only line ends the records, whether or not it is empty
    for tail in ("\n \n", "\n \n\n", "\n\t\r\n \n", "\n  ", "\n\r", "\n\r\r\n", "\n\n \n")
        @test length(PAFReader(collect, IOBuffer(good * tail))) == 1
        reader = PAFReader(IOBuffer(good * tail))
        @test try_next!(reader) isa PAFRecord
        @test isempty(reader)
        @test try_next!(reader) === nothing
    end
    @test isempty(PAFReader(IOBuffer(" \n\t\n")))

    # ... after which only whitespace may follow
    for tail in ("\n \n" * good, "\n\t\nx", "\n \r\n\n" * good * "\n")
        reader = PAFReader(IOBuffer(good * tail))
        @test try_next!(reader) isa PAFRecord
        @test !isempty(reader)
        @test try_next!(reader).kind == Errors.EmptyNonTailingLine
    end

    # Trailing whitespace on a record is not stripped
    reader = PAFReader(IOBuffer(good * " \n"))
    @test try_next!(reader).kind == Errors.BadInteger
    reader = PAFReader(IOBuffer(good * "\tab:Z:x \n" * good * "\t\n"))
    @test aux_data(try_next!(reader))["ab"] == "x "
    @test try_next!(reader) isa PAFRecord # empty aux field is only invalid when accessed

    # A line with non-whitespace after leading whitespace is an ordinary error
    reader = PAFReader(IOBuffer(good * "\n \tx\n"))
    @test try_next!(reader) isa PAFRecord
    @test try_next!(reader).kind == Errors.BadInteger
end

@testset "Skipping invalid lines" begin
    good = "my_qname\t301156\t3\t301145\t+\tmy_tname\t6701780\t2764860\t3066002\t301142\t301142\t0"
    kinds(results) = [r isa PAFRecord ? :record : (r.kind, r.line) for r in results]
    function read_all(reader)
        results = []
        while (res = try_next!(reader)) !== nothing
            push!(results, res)
        end
        return results
    end

    # try_next! advances past invalid lines
    data = join([good, "bad", good, good * "\tx", good * " ", good], '\n')
    @test kinds(read_all(PAFReader(IOBuffer(data)))) ==
        [:record, (Errors.TooFewFields, 2), :record, :record, (Errors.BadInteger, 5), :record]

    # Including a truncated last line, as from a killed job
    reader = PAFReader(IOBuffer(good * "\n" * good[1:30]))
    @test try_next!(reader) isa PAFRecord
    @test try_next!(reader).line == 2
    @test try_next!(reader) === nothing
    @test isempty(reader)

    # A run of whitespace-only lines followed by non-whitespace is one error,
    # reported at the first whitespace-only line. Reading continues after the run.
    data = good * "\n\n \r\n\t\n" * good * "\n\nbad\n\n"
    @test kinds(read_all(PAFReader(IOBuffer(data)))) == [
        :record, (Errors.EmptyNonTailingLine, 2), :record,
        (Errors.EmptyNonTailingLine, 6), (Errors.TooFewFields, 7),
    ]

    # Iteration throws, but also advances, so it can be resumed after catching
    reader = PAFReader(IOBuffer(join([good, "bad", good], '\n')))
    @test iterate(reader) isa Tuple
    @test_throws PairwiseMappingFormat.ParserException iterate(reader)
    @test iterate(reader) isa Tuple
    @test iterate(reader) === nothing

    # Mixing iteration and try_next! gives the same results
    reader = PAFReader(IOBuffer(join(["bad", "bad", good], '\n')))
    @test_throws PairwiseMappingFormat.ParserException iterate(reader)
    @test try_next!(reader).line == 2
    @test only(collect(reader)).qname == "my_qname"
end

@testset "isempty" begin
    good = "my_qname\t301156\t3\t301145\t+\tmy_tname\t6701780\t2764860\t3066002\t301142\t301142\t0"

    # isempty does not consume records
    reader = PAFReader(IOBuffer(join(fill(good, 3), '\n')))
    @test !isempty(reader)
    @test !isempty(reader)
    @test length(collect(reader)) == 3
    @test isempty(reader)

    # isempty does not invalidate records iterated with copy=true, even when it needs to
    # fill the buffer
    expected = parse(PAFRecord, good)
    for chunks in ([good * "\n", good * "\n"], [good * "\n", "\n", "\r\n"], [good])
        reader = PAFReader(ChunkReader(chunks))
        (rec, _) = iterate(reader)
        while !isempty(reader)
            @test rec.qname == expected.qname
            @test rec.tname == expected.tname
            @test aux_data(rec) == aux_data(expected)
            (rec, _) = iterate(reader)
        end
        @test rec.qname == expected.qname
        @test rec.tname == expected.tname
        @test iterate(reader) === nothing
    end

    # isdone does not fill the buffer when a view returned by try_next! is live,
    # but isempty does
    reader = PAFReader(ChunkReader([good * "\n", good * "\n"]); copy = false)
    rec = try_next!(reader)
    @test Base.isdone(reader) === missing
    @test rec.qname == "my_qname"
    @test !isempty(reader)
    @test try_next!(reader) isa PAFRecord
    @test isempty(reader)

    inputs = [
        "", "\n", "\r\n", "\n\n \t\r\n", good, good * "\n", good * "\r\n\n",
        good * "\n\n" * good, good * "\nbad", good * "\n \n", good * "\n\r",
        good * "\n\r\r\n", " ", "\nV", good * "\n\n\n\nx",
        good * "\n \n" * good, good * "\n\t \r\n  ", " \tx", good * "\n \n\tx",
        "bad\n" * good, good * "\nbad\n" * good * "\n", "bad\r\nbad\n\n \n" * good,
    ]
    # Empty chunks are not allowed, since fill_buffer returning zero signals EOF
    function chunkings(data)
        splits = [filter!(!isempty, [data[1:i], data[(i + 1):end]]) for i in 0:length(data)]
        return push!(splits, [string(c) for c in data])
    end
    # Summarize the results of reading all of `reader`, continuing past errors.
    # Each call must advance the reader, so the number of results is bounded.
    function read_all(reader, max_results)
        results = []
        while (res = try_next!(reader)) !== nothing
            push!(results, res isa PAFRecord ? copy(getfield(res, :data)) : res)
            length(results) > max_results && error("Reader did not advance")
        end
        return results
    end
    for data in inputs
        expected_results = read_all(PAFReader(IOBuffer(data)), ncodeunits(data))
        for chunks in chunkings(data)
            @test read_all(PAFReader(ChunkReader(chunks)), ncodeunits(data)) == expected_results
        end
    end

    for data in inputs, chunks in chunkings(data), copy in (true, false)
        # isdone and isempty agree with try_next!, and isdone never invalidates
        # the record returned by try_next!. With copy=true, no view is live, so
        # isdone is never missing, and isempty does not invalidate the record either.
        # Reading continues past errors.
        reader = PAFReader(ChunkReader(chunks); copy)
        (rec, expected) = (nothing, nothing)
        while true
            maybe_done = Base.isdone(reader)
            copy && @test maybe_done isa Bool
            rec === nothing || @test getfield(rec, :data) == getfield(expected, :data)
            done = isempty(reader)
            copy && rec !== nothing && @test getfield(rec, :data) == getfield(expected, :data)
            maybe_done === missing || @test maybe_done == done
            res = try_next!(reader)
            @test done == (res === nothing)
            res === nothing && break
            (rec, expected) = res isa PAFRecord ? (res, Base.copy(res)) : (nothing, nothing)
        end

        # When iterating with copy=true, no view is live, so isdone is never missing
        reader = PAFReader(ChunkReader(chunks))
        while true
            done = Base.isdone(reader)
            @test done isa Bool
            it = try
                iterate(reader)
            catch e
                e isa PairwiseMappingFormat.ParserException || rethrow()
                e
            end
            @test done == (it === nothing)
            it === nothing && break
        end
    end
end

@testset "Aqua" begin
    Aqua.test_all(PairwiseMappingFormat; ambiguities = false)
end

end # module
