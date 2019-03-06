using LinearAlgebra, SparseArrays, InfiniteArrays, FillArrays, LazyArrays, Statistics, DSP, Test
    import InfiniteArrays: OrientedInfinity, OneToInf, InfUnitRange, InfStepRange


@testset "∞" begin
    @test ∞ ≠ 1
    @test ∞ == ∞
    @test ∞ == Inf

    @test isless(1, ∞)
    @test !isless(Inf, ∞)
    @test !isless(∞, Inf)
    @test !isless(∞, 1)

    @test !isless(∞, ∞)
    @test !(∞ < ∞)
    @test ∞ ≤ ∞
    @test !(∞ > ∞)
    @test ∞ ≥ ∞

    @test ∞ + ∞ ≡ ∞
    @test ∞ + 1 ≡ ∞
    @test *(∞) ≡ ∞
    @test ∞*∞ ≡ ∞

    # oriented infinity
    @test OrientedInfinity(∞) ≡ convert(OrientedInfinity, ∞) ≡ OrientedInfinity() ≡
        OrientedInfinity(false)

    @test -∞ ≡ OrientedInfinity(true)
    @test +∞ ≡ ∞

    @test ∞ == +∞ == OrientedInfinity(∞)
    @test ∞ ≠  -∞
    @test 1-∞ == -∞

    @test (-∞)*(-∞) ≡ ∞*OrientedInfinity(∞) ≡ OrientedInfinity(∞)*∞

    @test  isless(-∞, 1)
    @test !isless(-∞, -Inf)
    @test !isless(-Inf, -∞)
    @test !isless(1, -∞)

    @test max(1,∞) == max(∞,1) == ∞
    @test min(1,∞) == min(∞,1) == 1
    @test maximum([1,∞]) == ∞
    @test minimum([1,∞]) == 1

    @test OrientedInfinity(true) + OrientedInfinity(true) == OrientedInfinity(true)
    @test OrientedInfinity(false) + OrientedInfinity(false) == OrientedInfinity(false)
    @test OrientedInfinity(true)+1 == OrientedInfinity(true)
    @test OrientedInfinity(false)+1 == OrientedInfinity(false)



    @test exp(im*π/4)*∞ == Inf+im*Inf
    @test exp(im*π/4)+∞ == ∞
end

@testset "ranges" begin
    @test size(10:1:∞) == (∞,)
    @testset "colon" begin
        @test @inferred(10:1:∞) === @inferred(range(10; step=1, length=∞))
        @inferred(1:0.2:∞)  === @inferred(range(1; step=0.2, length=∞))
        @inferred(1.0:0.2:∞)  === @inferred(range(1.0; step=0.2, length=∞))
        @inferred(1:∞) === @inferred(range(1; length=∞))
    end
    @test_throws ArgumentError 2:-.2:∞
    @test_throws ArgumentError 0.0:-∞

    @testset "indexing" begin
        L32 = @inferred(Int32(1):∞)
        L64 = @inferred(Int64(1):∞)
        @test @inferred(L32[1]) === Int32(1) && @inferred(L64[1]) === Int64(1)
        @test L32[2] == 2 && L64[2] == 2
        @test L32[3] == 3 && L64[3] == 3
        @test L32[4] == 4 && L64[4] == 4

        @test @inferred((1.0:∞)[1]) === 1.0
        @test @inferred((1.0f0:∞)[1]) === 1.0f0
        @test @inferred((Float16(1.0):∞)[1]) === Float16(1.0)

        @test @inferred((0.1:0.1:∞)[2]) === 0.2
        @test @inferred((0.1f0:0.1f0:∞)[2]) === 0.2f0
        @test @inferred((1:∞)[1:4]) === 1:4
        @test @inferred((1.0:∞)[1:4]) === 1.0:4
        @test (2:∞)[1:4] == 2:5
        @test (1:∞)[2:5] === 2:5
        @test (1:∞)[2:2:5] === 2:2:4
        @test (1:2:∞)[2:6] === 3:2:11
        @test (1:2:∞)[2:3:7] === 3:6:13

        @test isempty((1:∞)[5:4])
        @test_throws BoundsError (1:∞)[8:-1:-2]
    end
    @testset "length" begin
        @test length(.1:.1:∞) == ∞
        @test length(1.1:1.1:∞) == ∞
        @test length(1.1:1.3:∞) == ∞
        @test length(1:1:∞) == ∞
        @test length(1:.2:∞) == ∞
        @test length(1.:.2:∞) == ∞
        @test length(2:-.2:-∞) == ∞
        @test length(2.:-.2:-∞) == ∞
    end

    @test_throws ArgumentError (2:.2:-∞)
    @test_throws ArgumentError (2.:.2:-∞)
    @test_throws ArgumentError (1:-∞)


    @testset "intersect" begin
        @test intersect(1:∞, 2:3) == 2:3
        @test intersect(2:3, 1:∞) == 2:3
        @test intersect(1:∞, 2:∞) == 2:∞
        @test intersect(2:∞, 1:∞) == 2:∞
        @test intersect(-3:∞, 2:8) == 2:8
        @test intersect(1:∞, -2:3:15) == 1:3:15
        @test intersect(1:∞, -2:3:∞) == 1:3:∞
        @test intersect(1:11, -2:2:∞) == intersect(-2:2:∞,1:11) == 2:2:10
        @test intersect(1:∞, -2:1:15) == 1:15
        @test intersect(1:∞, 15:-1:-2) == 1:15
        @test intersect(1:∞, 15:-4:-2) == 3:4:15
        @test intersect(-20:∞, -10:3:-2) == -10:3:-2
        @test isempty(intersect(-5:5, -6:13:∞))
        @test isempty(intersect(1:∞, 15:4:-2))

        @test @inferred(intersect(0:3:∞, 0:4:∞)) == intersect(0:4:∞, 0:3:∞) == 0:12:∞

        @test intersect(24:-3:0, 0:4:∞) == 0:12:24
        @test_throws ArgumentError intersect(1:6:∞, 0:4:∞) # supporting empty would break type inferrence

        @test intersect(1:∞,3) == 3:3
        @test intersect(1:∞, 2:∞, UnitRange(3,7), UnitRange(4,6)) == UnitRange(4,6)

        @test intersect(1:∞, 2) === intersect(2, 1:∞) === 2:2
        @test_skip intersect(1.0:∞, 2) == intersect(2, 1.0:∞) == [2.0] # gives infinite loop
    end

    @testset "sort/sort!/partialsort" begin
        @test sort(1:∞) == 1:∞
        @test sort!(1:∞) == 1:∞
    end
    @testset "in" begin
        @test 0 in UInt(0):100:∞
        @test 0 in 0:-100:-∞
        @test ∞ ∉ UInt(0):100:∞

        @test !(3.5 in 1:∞)
        @test (3 in 1:∞)
        @test (3 in 5:-1:-∞)

        let r = 0.0:0.01:∞
            @test (r[30] in r)
        end
        let r = (-4*Int64(maxintfloat(Int === Int32 ? Float32 : Float64))):∞
            @test (3 in r)
            @test (3.0 in r)
        end
    end
    @testset "in() works across types, including non-numeric types (#21728)" begin
        @test 1//1 in 1:∞
        @test 1//1 in 1.0:∞
        @test !(5//1 in 6:∞)
        @test !(5//1 in 6.0:∞)
        @test Complex(1, 0) in 1:∞
        @test Complex(1, 0) in 1.0:∞
        @test Complex(1.0, 0.0) in 1:∞
        @test Complex(1.0, 0.0) in 1.0:∞
        @test_skip !(Complex(1, 1) in 1:∞)  # this is an infinite-loop at the moment
        @test_skip !(Complex(1, 1) in 1.0:∞)
        @test_skip !(Complex(1.0, 1.0) in 1:∞)
        @test_skip !(Complex(1.0, 1.0) in 1.0:∞)
        @test !(π in 1:∞)
    end


    @testset "indexing range with empty range (#4309)" begin
        @test (3:∞)[5:4] == 7:6
        @test (0:2:∞)[7:6] == 12:2:10
    end

    @testset "indexing with negative ranges (#8351)" begin
        for a=[3:∞, 0:2:∞], b=[0:1, 2:-1:0]
            @test_throws BoundsError a[b]
        end
    end

    @testset "sums of ranges" begin
        @test sum(1:∞) ≡ mean(1:∞) ≡ median(1:∞) ≡ ∞
        @test sum(0:∞) ≡ mean(1:∞) ≡ median(1:∞) ≡ ∞
        @test sum(0:2:∞) ≡ mean(0:2:∞) ≡ median(0:2:∞) ≡ OrientedInfinity(∞)
        @test sum(0:-2:-∞) ≡ mean(0:-2:-∞) ≡ median(0:-2:-∞) ≡ -∞
    end

    @testset "broadcasted operations with scalars" begin
        @test broadcast(-, 1:∞, 2) ≡ -1:∞
        @test broadcast(-, 1:∞, 0.25) ≡ 1-0.25:∞
        @test broadcast(+, 1:∞, 2) ≡ 3:∞
        @test broadcast(+, 1:∞, 0.25) ≡ 1+0.25:∞
        @test broadcast(+, 1:2:∞, 1) ≡ 2:2:∞
        @test broadcast(+, 1:2:∞, 0.3) ≡ 1+0.3:2:∞
        @test broadcast(-, 1:2:∞, 1) ≡ 0:2:∞
        @test broadcast(-, 1:2:∞, 0.3) ≡ 1-0.3:2:∞
        @test broadcast(-, 2, 1:∞) ≡ 1:-1:-∞
    end

    @testset "near-equal ranges" begin
        @test 0.0:0.1:∞ != 0.0f0:0.1f0:∞
    end

    @testset "comparing InfiniteUnitRanges and OneToInf" begin
        @test 1:2:∞ == 1:2:∞ != 1:3:∞ != 2:3:∞ == 2:3:∞ != 2:∞
        @test 1:1:∞ == 1:∞ == 1:∞ == OneToInf() == OneToInf()
    end

    @testset "Base.OneTo (misleading) overrides" begin
        @test Base.OneTo{BigInt}(∞) isa OneToInf{BigInt}
        @test Base.OneTo(∞) isa OneToInf{Int}
    end

    @testset "issue #6973" begin
        r1 = 1.0:0.1:∞
        r2 = 1.0f0:0.2f0:∞
        r3 = 1:2:∞
        @test r1 + r1 == 2*r1
        @test_broken r1 + r2 == 2.0:0.3:∞
        @test (2r1)-3r1 == -1:(2step(r1)-3step(r1)):-∞
        @test_broken (r1 + r2) - r2 == r1
        @test r1 + r3 == 2.0:2.1:∞
        @test r3 + r3 == 2 * r3
    end

    @testset "Preservation of high precision upon addition" begin
        let r = (-0.1:0.1:∞) + broadcast(+, -0.3:0.1:∞, 1e-12)
            @test_broken r[3] == 1e-12
        end
    end

    @testset "issue #8584" begin
        @test (0:1//2:∞)[1:2:3] == 0:1//1:1
    end

    @testset "issue #9962" begin
        @test eltype(0:1//3:∞) <: Rational
        @test (0:1//3:∞)[1] == 0
        @test (0:1//3:∞)[2] == 1//3
    end
    @testset "converting ranges (issue #10965)" begin
        @test promote(0:∞, UInt8(2):∞) === (0:∞, 2:∞)
        @test convert(InfUnitRange{Int}, 0:∞) === 0:∞
        @test convert(InfUnitRange{Int128}, 0:∞) === Int128(0):∞

        @test InfUnitRange{Int16}(1:∞) ≡ AbstractVector{Int16}(1:∞) ≡
                AbstractArray{Int16}(1:∞) ≡ Int16(1):∞

        @test OneToInf{Int16}(OneToInf()) ≡ AbstractVector{Int16}(OneToInf()) ≡
                AbstractArray{Int16}(OneToInf()) ≡ OneToInf{Int16}()

        @test promote(0:1:∞, UInt8(2):UInt8(1):∞) === (0:1:∞, 2:1:∞)
        @test convert(InfStepRange{Int,Int}, 0:1:∞) === 0:1:∞
        @test convert(InfStepRange{Int128,Int128}, 0:1:∞) === Int128(0):Int128(1):∞

        @test promote(0:1:∞, 2:∞) === (0:1:∞, 2:1:∞)
        @test convert(InfStepRange{Int128,Int128}, 0:∞) === Int128(0):Int128(1):∞
        @test convert(InfStepRange, 0:∞) === 0:1:∞
        @test convert(InfStepRange{Int128,Int128}, 0.:∞) === Int128(0):Int128(1):∞

        @test_broken promote(0f0:inv(3f0):∞, 0.:2.:∞) === (0:1/3:∞, 0.:2.:∞)

        @test promote(0:1/3:∞, 0:∞) === (0:1/3:∞, 0.:1.:∞)
    end

    @testset "inf-range[inf-range]" begin
        @test (1:∞)[1:∞] == 1:∞
        @test (1:∞)[OneToInf()] == 1:∞
        @test (1:∞)[2:∞] == 2:∞
        @test_throws BoundsError (1:∞)[-1:∞]
        @test (1:-1:-∞)[1:∞] == 1:-1:-∞
    end

    @testset "OneToInf" begin
        let r = OneToInf()
            @test !isempty(r)
            @test length(r) == ∞
            @test size(r) == (∞,)
            @test step(r) == 1
            @test first(r) == 1
            @test last(r) == ∞
            @test minimum(r) == 1
            @test maximum(r) == ∞
            @test r[2] == 2
            @test r[2:3] === 2:3
            @test_throws BoundsError r[0]
            @test broadcast(+, r, 1) === 2:∞
            @test 2*r === 2:2:∞
            @test r + r === 2:2:∞

            @test r - r === Zeros{Int}(∞)

            @test intersect(r, Base.OneTo(2)) == Base.OneTo(2)
            @test intersect(r, 0:5) == 1:5
            @test intersect(r, 2) === intersect(2, r) === 2:2
        end
    end

    @testset "show" begin
        # NOTE: Interpolating Int to ensure it's displayed properly across 32- and 64-bit
        @test summary(1:∞) == "InfUnitRange{$Int} with indices OneToInf()"
        @test Base.inds2string(axes(1:∞)) == "OneToInf()"
    end

    @testset "end" begin
        @test Base.OneTo(∞)[end] ≡ Base.OneTo(∞)[∞] ≡ ∞
        @test (1:∞)[end] ≡ (1:∞)[∞] ≡ ∞
        @test (1:2:∞)[end] ≡ (1:2:∞)[∞] ≡ ∞
        @test (1.0:2:∞)[end] ≡ (1.0:2:∞)[∞] ≡ ∞
    end
end

@testset "fill" begin
    for A in (Zeros(∞), Fill(1,∞), Ones(∞))
        @test length(A) == ∞
    end
    @test size(Zeros(∞,5)) === (∞,5)
    @test size(Zeros(5,∞)) === (5,∞)
end

@testset "diagonal" begin
    D = Diagonal(1:∞)
    @test D[1:10,1:10] == Diagonal(1:10)
end


@testset "concat" begin
    A = Vcat(1:10, 1:∞)
    @test @inferred(length(A)) == ∞
    @test @inferred(A[5]) == A[15] == 5
    @test A[end] == @inferred(A[∞]) == ∞
    @test_throws BoundsError Vcat(1:10)[∞]

    A = Vcat(Ones(1,∞), Zeros(2,∞))
    @test @inferred(size(A)) == (3,∞)
    @test @inferred(A[1,5]) == 1
    @test @inferred(A[3,5]) == 0
    @test_throws BoundsError A[4,1]

    A = Vcat(Ones{Int}(1,∞), Diagonal(1:∞))
    @test @inferred(size(A)) ≡ (∞,∞)
    @test @inferred(A[1,5]) ≡ 1
    @test @inferred(A[5,5]) ≡ 0
    @test @inferred(A[6,5]) ≡ 5
    @test_throws BoundsError A[-1,1]

    A = Vcat(Ones{Float64}(1,∞), Diagonal(1:∞))
    @test @inferred(size(A)) ≡ (∞,∞)
    @test @inferred(A[1,5]) ≡ 1.0
    @test @inferred(A[5,5]) ≡ 0.0
    @test @inferred(A[6,5]) ≡ 5.0
    @test_throws BoundsError A[-1,1]

    A = Vcat(1, Zeros(∞))
    @test @inferred(A[1]) ≡ 1.0
    @test @inferred(A[2]) ≡ 0.0

    A = Hcat(Ones(∞), Zeros(∞,2))
    @test @inferred(size(A)) == (∞,3)
    @test @inferred(A[5,1]) == 1
    @test @inferred(A[5,3]) == 0
    @test_throws BoundsError A[1,4]

    A = Hcat(Ones{Int}(∞), Diagonal(1:∞))
    @test @inferred(size(A)) ≡ (∞,∞)
    @test @inferred(A[5,1]) ≡ 1
    @test @inferred(A[5,5]) ≡ 0
    @test @inferred(A[5,6]) ≡ 5
    @test_throws BoundsError A[-1,1]
end

# This should be generalized, but it at the moment
# it is restricted to a single Number. Support smart
# addition for any number of Number/SVector's would be better
# allowibng for the tail to be variable lenth
@testset "Vcat special case" begin
    @test Vcat(1,Zeros{Int}(∞)) + Vcat(3,Zeros{Int}(∞)) ≡
          Vcat(1,Zeros{Int}(∞)) .+ Vcat(3,Zeros{Int}(∞)) ≡
          Vcat(4,Zeros{Int}(∞))

          size(Vcat(1:∞)) ≡ (∞,)
end

@testset "Fill indexing" begin
    B = Ones(∞,∞)
    @test IndexStyle(B) == IndexCartesian()
    V = view(B,:,1)
    @test_broken size(V) == (∞,1)
    V = view(B,1,:)
    @test size(V) == (∞,)
    V = view(B,1:1,:)
    @test size(V) == (1,∞)
end

@testset "∞ BroadcastArray" begin
    A = 1:∞
    B = BroadcastArray(exp, A)
    @test length(B) == ∞
    @test B[6] == exp(6)
    @test exp.(A) ≡ B
    B = Diagonal(1:∞) .+ 1
    @test B isa BroadcastArray{Int}
    @test B[1,5] ≡ 1
    @test B[6,6] == 6+1
    B = Diagonal(1:∞) - Ones{Int}(∞,∞) # lowers to broadcast
    @test B isa BroadcastArray{Int}
    @test B[1,5] ≡ -1
    @test B[6,6] == 6-1
end

@testset "Taylor ODE" begin
    e₁ = Vcat(1, Zeros(∞))
    D = Hcat(Zeros(∞), Diagonal(1:∞))
    I_inf = Eye(∞)
    @test I_inf isa Eye{Float64,OneToInf{Int}}
    @test axes(I_inf) == (OneToInf{Int}(), OneToInf{Int}())
    @test eltype(I_inf) == Float64
    @test Base.BroadcastStyle(typeof(I_inf)) == LazyArrays.LazyArrayStyle{2}()
    L = Vcat(e₁', I_inf + D)
    @test L[1:3,1:3] == [1.0 0.0 0.0;
                         1.0 1.0 0.0;
                         0.0 1.0 2.0]
end

@testset "Cumsum and diff" begin
    @test cumsum(Ones(∞)) ≡ 1.0:1.0:∞
    @test cumsum(Fill(2,∞)) ≡ 2:2:∞
    @test cumsum(Ones{Int}(∞)) ≡ Base.OneTo(∞)
    @test cumsum(Ones{BigInt}(∞)) ≡ Base.OneTo{BigInt}(∞)

    @test diff(Base.OneTo(∞)) ≡ Ones{Int}(∞)
    @test diff(1:∞) ≡ Fill(1,∞)
    @test diff(1:2:∞) ≡ Fill(2,∞)
    @test diff(1:2.0:∞) ≡ Fill(2.0,∞)

    x = Vcat([3,4], Ones{Int}(5), 3, Fill(2,∞))
    y = @inferred(cumsum(x))
    @test y isa Vcat
    @test y[1:12] == cumsum(x[1:12])

    @test diff(x[1:10]) == diff(x)[1:9]
    @test diff(y)[1:20] == x[2:21]

    @test cumsum(x).arrays[2] ≡ 8:12
    @test last(y.arrays) == sum(x[1:9]):2:∞
    r = (3:4:∞)
    @test cumsum(r)[1:20] == cumsum(r[1:20])
end

@testset "Sub-array" begin
    @test Ones(∞)[3:∞] ≡ Ones(∞)
    @test Ones{Int}(∞)[4:6] ≡ Ones{Int}(3)
    @test (1:∞)[3:∞] ≡ 3:∞
end

@testset "Broadcast Fill Lowers" begin
    @test broadcast(+, Zeros{Int}(∞) , Fill(1,∞)) isa Fill
    @test broadcast(+, Zeros{Int}(∞) , Zeros(∞)) isa Zeros
    @test broadcast(*, Ones(∞), Ones(∞)) ≡ Ones(∞)
end

@testset "maximum/minimum Vcat" begin
    x = Vcat(1:2, [1,1,1,1,1], 3, Fill(4,∞))
    @test maximum(x) == 4
    @test minimum(x) == 1

    @test_throws ArgumentError maximum(exp.(1:∞))
end

@testset "conv" begin
    @test conv(1:∞, [2]) ≡ conv([2], 1:∞) ≡ 2:2:∞
    @test conv(1:2:∞, [2]) ≡ conv([2], 1:2:∞) ≡ 2:4:∞
    @test conv(1:∞, Ones(∞))[1:5] == conv(Ones(∞),1:∞)[1:5] == [1,3,6,10,15]
    @test conv(Ones(∞), Ones(∞)) ≡ 1.0:1.0:∞
    @test conv(Ones{Int}(∞), Ones{Int}(∞)) ≡ Base.OneTo(∞)
    @test conv(Ones{Bool}(∞), Ones{Bool}(∞)) ≡ Base.OneTo(∞)
end

@testset "show" begin
    @test repr(Vcat(1:∞)) == "[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, …]"
    @test repr(Vcat(2,1:∞)) == "[2, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, …]"
end

@testset "Array constructor errors" begin
    @test_throws ArgumentError Array{Float64}(undef, ∞)
    @test_throws ArgumentError Array{Float64}(undef, ∞, ∞)
    @test_throws ArgumentError Array{Float64}(undef, 1, ∞)
    @test_throws ArgumentError Array{Float64}(undef, ∞, 1)

    @test_throws ArgumentError Vector{Float64}(undef, ∞)
    @test_throws ArgumentError Matrix{Float64}(undef, ∞, ∞)
    @test_throws ArgumentError Matrix{Float64}(undef, 1, ∞)
    @test_throws ArgumentError Matrix{Float64}(undef, ∞, 1)
end
