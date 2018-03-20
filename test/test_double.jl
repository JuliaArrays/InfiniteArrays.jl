function test_doubly_infinite_arrays()
    @testset "doubly infinite sequences" begin
        @testset "zero padding" begin
            test_zero_padding()
        end
        @testset "periodic extension" begin
            test_periodic_extension()
        end
        @testset "symmetric extension" begin
            test_symmetric_extension()
        end
    end
end

function test_zero_padding()
    v1 = [1,2,3]
    # First, use a default offset and test basic functionality
    z1 = zero_padding(v1)
    @test eltype(z1) == eltype(v1)
    @test sum(abs.([z1[i]-v1[i] for i in eachindex(v1)])) == 0
    @test support(z1) == 1:3
    @test z1[0] == 0
    @test z1[4] == 0
    @test z1[2:3] == v1[2:3]

    z2 = reverse(z1)
    @test sum(abs.([z2[-i]-v1[i] for i in eachindex(v1)])) == 0

    # Test a different offset
    offset = 2
    z3 = CompactSequence(v1, offset)
    @test sum(abs.([z3[offset+i-1]-v1[i] for i in eachindex(v1)])) == 0
    supp = InfiniteArrays.data_support(z1, z3)
    @test first(supp) == 1
    @test last(supp) == length(v1)+1

    s = 4
    z4 = shift(z1, s)
    @test sum(abs.([z4[i+s]-v1[i] for i in eachindex(v1)])) == 0

    # Arithmetic
    z1 = CompactSequence(v1)
    z2 = CompactSequence(v1, 2)
    s1 = support(z1)
    a = z1+z2
    @test first(support(a)) == min(first(support(z1)), first(support(z2)))
    @test last(support(a)) == max(last(support(z1)), last(support(z2)))
    supp = support(a)
    @test maximum(abs.(a[supp] - z1[supp] - z2[supp])) == 0
    b = z1-z2
    @test maximum(abs.(b[supp] - z1[supp] + z2[supp])) == 0
    c = 2*z1
    @test maximum(abs.(c[s1]-2*z1[s1])) == 0
    d = -z1
    @test maximum(abs.(d[s1]+z1[s1])) == 0
    e = z1/2
    @test maximum(abs.(e[s1]-z1[s1]/2)) <= 10eps(eltype(e))
end

function test_periodic_extension()
    v1 = [1,2,3,4,5]
    z1 = PeriodicExtension(v1)
    @test eltype(z1) == eltype(v1)
    @test sum(abs.([z1[i]-v1[i] for i in eachindex(v1)])) == 0
    @test z1[0] == v1[end]
    @test z1[6] == v1[1]
    @test z1[2:3] == v1[2:3]

    # Arithmetic
    v2 = [3,4,2,1,6]
    z1 = PeriodicExtension(v1)
    z2 = PeriodicExtension(v2, 2)
    p = period(z1)
    a = z1+z2
    @test period(a) == period(z1) == period(z2)
    s = eachindex(v2)
    @test maximum(abs.(a[s] - z1[s] - z2[s])) == 0
    b = z1-z2
    @test maximum(abs.(b[s] - z1[s] + z2[s])) == 0
    c = 2*z1
    @test maximum(abs.(c[s]-2*z1[s])) == 0
    d = -z1
    @test maximum(abs.(d[s]+z1[s])) == 0
    e = z1/2
    @test maximum(abs.(e[s]-z1[s]/2)) <= 10eps(eltype(e))
end

function test_symmetric_extension()
end
