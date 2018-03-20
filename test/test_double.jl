function test_doubly_infinite_arrays()
    @testset "doubly infinite sequences" begin
        @testset "zero padding" begin
            test_zero_padding()
        end
    end
end

function test_zero_padding()
    v1 = [1,2,3]
    # First, use a default offset and test basic functionality
    z1 = zero_padding(v1)
    @test eltype(z1) == eltype(v1)
    @test sum(abs.([z1[i]-v1[i] for i in eachindex(v1)])) == 0
    @test leftmost(z1) == 1
    @test rightmost(z1) == 3
    @test z1[0] == 0
    @test z1[4] == 0
    @test z1[2:3] == v1[2:3]

    z2 = reverse(z1)
    @test sum(abs.([z2[-i]-v1[i] for i in eachindex(v1)])) == 0

    # Test a different offset
    offset = 2
    z3 = CompactSequence(v1, offset)
    @test sum(abs.([z3[offset+i-1]-v1[i] for i in eachindex(v1)])) == 0
end
