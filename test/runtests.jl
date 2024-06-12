module TestPAF

using PairwiseMappingFormat
using PairwiseMappingFormat: Errors, try_parse
using Test
using FormatSpecimens
using MemViews

using SAMAuxData.SAM: AuxTag

function cmp_aux(record::PAFRecord, other)
    aux = sort!(collect(aux_data(record)); by=first)
    other = sort!(collect(other); by=first)
    @test length(aux) == length(other)
    for ((kaux, vaux), (ko, vo)) in zip(aux, other)
        @test kaux == AuxTag(ko)
        if vaux isa Union{Real, AbstractString, Char}
            @test vaux == vo
        elseif vaux isa AbstractVector
            @test eltype(vaux) == eltype(vo)
            @test vaux == vo
        else
            @warn typeof(vaux)
            @test false
        end
    end
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
	end

	@testset "Aux data" begin
	    cmp_aux(rec1, ["tp" => 'P', "cm" => 29990, "s1" => 301142, "s2" => 301142, "dv" => 0, "rl" => 0])
	    cmp_aux(rec2, ["tp" => 'P', "cm" => 29957, "s1" => 299269, "s2" => 299269, "dv" => 0, "rl" => 38])
	    cmp_aux(rec3, ["tp" => 'P', "cm" => 28871, "s1" => 288618, "s2" => 288618, "dv" => 0, "rl" => 57])
    end
end

@testset "Edge case records" begin
    (rec1, rec2, rec3) = PAFReader(collect, open(joinpath(path_of_format("PAF"), "good2.paf")))
    @test rec1.qname == "中国科学院大学"
    @test rec1.qlen == typemax(Int32)
    cmp_aux(rec1, ["tp" => 'S', "cm" => 29990, "s1" => 301142, "dv" => -13.2442f-18, "rl" => 0, "xX" => [0xea, 0x19, 0xfe, 0x1b], "K1" => Int8[15, 9, 22, 127, -12, -128, 0]])
    cmp_aux(rec2, ["dv" => 1.243f-4, "rl" => 38, "kk" => "!=[]! And then~~  "])
    cmp_aux(rec3, ["dv" => +17f21, "rl" => 38])
end

function with_replaced(s::String, field::Integer, new::String)
    fields = split(s, '\t')
    fields[field] = new
    y = try_parse(join(fields, "\t"))
    y isa PairwiseMappingFormat.ParserError ? y.kind : y
end

@testset "Errors" begin
    good = "my_qname	301156	3	301145	+	my_tname	6701780	2764860	3066002	301142	301142	0"
    @test with_replaced(good, 1, "my other name") isa PAFRecord
    
    # Negative numbers, zero numbers
    @test with_replaced(good, 2, "-1") == Errors.BadInteger

    for allow_zero in [3, 8, 10, 12]
        @test with_replaced(good, allow_zero, "0") isa PAFRecord
    end

    for not_allow_zero in [2, 4, 7, 9, 11]
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
        @test y isa PairwiseMappingFormat.ParserError && y.kind == Errors.TooFewFields
    end

    # Missing field
    for int_field in [2, 3, 4, 7, 8, 9, 10, 11, 12]
        @test with_replaced(good, int_field, "") == Errors.EmptyInteger
    end
    
    # Mapq too high
    @test with_replaced(good, 12, "256") == Errors.IntegerOverflow
end

# We assume that aux data is correctly tested in its own package, here we just test
# that the aux data has the right span and byte content.
function test_aux_content(record_data::String, target_aux::String)
    rec = parse(PAFRecord, record_data)
    @test collect(MemView(aux_data(rec))) == codeunits(target_aux)
end

@testset "Auxiliary data" begin
    good = "my_qname	301156	3	301145	+	my_tname	6701780	2764860	3066002	301142	301142	0	tp:A:P	kL:Z:325 1	s1:i:301142	x3:H:a40140	dv:f:0.0000	rl:i:0"
    goodaux = "tp:A:P	kL:Z:325 1	s1:i:301142	x3:H:a40140	dv:f:0.0000	rl:i:0"
    test_aux_content(good, goodaux)

    stringaux_with_whitespace = "\tHA:Z:lsjd   " 
    test_aux_content(good * stringaux_with_whitespace, goodaux * stringaux_with_whitespace)
end

# Note: We currently DONT trim whitespace off ends of records, for the sole reason that
# whitespace is allowed in string AuxData.
# This is very annoying, because it makes parsing less robust to erroneous trailing whitespace.
# However, I can't think of a good way to deal with this.
# TODO: Find a better method
@testset "Trailing whitespace" begin
    data = "my_qname	301156	3	301145	+	my_tname	6701780	2764860	3066002	301142	301142	0	an:Z: dljdl \r\n"
    data *= "some	5	3	4	+	other	9	5	8	2	2	0\r\n"

    reader = PAFReader(IOBuffer(data))
    (rec1, rec2) = collect(reader)

    @test aux_data(rec1)["an"] == " dljdl "
end

end # module
