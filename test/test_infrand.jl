@testset "InfRandVector" begin
    @testset "Default constructor" begin
        Random.seed!(123)
        seq = InfiniteArrays.InfRandVector()
        val = seq[1]
        @test seq[1] == val # constant 
        @test seq[1:10000] == seq[1:10000]
        @test seq[1] == val # didn't change after resizing
        @inferred seq[1]
        Random.seed!(123)
        _seq = [rand() for _ in 1:1000]
        @test seq[1:1000] == _seq[1:1000]
        @test size(seq) == (ℵ₀,)
        @test length(seq) == ℵ₀
        @test axes(seq) == (1:ℵ₀,)
    end

    @testset "Providing an RNG and a distribution" begin
        rng = MersenneTwister(123)
        seq = InfiniteArrays.InfRandVector(rng, Float16)
        rng2 = MersenneTwister(123)
        @test seq[1:10000] == [rand(rng2, Float16) for _ in 1:10000]
    end

    @testset "Distributions.jl" begin
        dist = Normal(0.3, 1.7) # do Normal{Float32} for example if you want that number type
        rng = Xoshiro(5)
        seq = InfiniteArrays.InfRandVector(rng, dist)
        rng2 = Xoshiro(5)
        @test seq[1:100] == [0.3 + 1.7randn(rng2) for _ in 1:100]
    end
end

@testset "InfRandMatrix" begin
    @testset "Default constructor" begin
        Random.seed!(123)
        A = InfiniteArrays.InfRandMatrix(5)
        @test size(A) == (5, ∞)
        @test axes(A) == (1:5, 1:∞)
        val = A[1, 1]
        @inferred A[5, 5]
        @test A[1:3, 1:100] == A[1:3, 1:100]
        @test_throws BoundsError A[0, 1]
        Random.seed!(123)
        _A = [rand() for _ in 1:5, _ in 1:1000]
        @test _A == A[1:5, 1:1000]
        @test (A+A)[1:5, 1:1000] ≈ 2_A
        @test (A'+A')[1:1000, 1:5] ≈ 2_A'
        @test A[11] ≈ A.seq[11] ≈ A[1, 3]
    end

    @testset "Providing an RNG and a distribution" begin
        rng = MersenneTwister(123)
        seq = InfiniteArrays.InfRandMatrix(rng, 10; dist=Float16)
        rng2 = MersenneTwister(123)
        @test seq[1:10000] == [rand(rng2, Float16) for _ in 1:10000]
    end

    @testset "Distributions.jl" begin
        dist = Normal(0.3, 1.7)
        rng = Xoshiro(5)
        seq = InfiniteArrays.InfRandMatrix(rng, 2; dist)
        rng2 = Xoshiro(5)
        @test seq[1:1000] == [0.3 + 1.7randn(rng2) for _ in 1:1000]
    end
end