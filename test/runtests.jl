using LinearAlgebra, SparseArrays, InfiniteArrays, Infinities, FillArrays, LazyArrays, Statistics, Test, Base64
using BandedMatrices
import InfiniteArrays: InfUnitRange, InfStepRange, OneToInf, NotANumber, oneto, unitrange
import LazyArrays: CachedArray, MemoryLayout, LazyLayout, DiagonalLayout, LazyArrayStyle, colsupport, DualLayout
import BandedMatrices: _BandedMatrix, BandedColumns
import Base.Broadcast: broadcasted, Broadcasted, instantiate

using Aqua
@testset "Project quality" begin
    Aqua.test_all(InfiniteArrays, ambiguities=false, piracies=false)
end

@testset "construction" begin
    @testset "Array constructor errors" begin
        @test_throws ArgumentError Array{Float64}(undef, ∞)
        @test_throws ArgumentError Array{Float64}(undef, ∞, ∞)
        @test_throws ArgumentError Array{Float64}(undef, 1, ∞)
        @test_throws ArgumentError Array{Float64}(undef, ∞, 1)

        @test_throws ArgumentError Vector{Float64}(undef, ∞)
        @test_throws ArgumentError Matrix{Float64}(undef, ∞, ∞)
        @test_throws ArgumentError Matrix{Float64}(undef, 1, ∞)
        @test_throws ArgumentError Matrix{Float64}(undef, ∞, 1)

        @test_throws ArgumentError Array{Float64}(undef, (∞,))
        @test_throws ArgumentError Array{Float64}(undef, (∞, ∞))
        @test_throws ArgumentError Array{Float64}(undef, (1, ∞))
        @test_throws ArgumentError Array{Float64}(undef, (∞, 1))

        @test_throws ArgumentError Vector{Float64}(undef, (∞,))
        @test_throws ArgumentError Matrix{Float64}(undef, (∞, ∞))
        @test_throws ArgumentError Matrix{Float64}(undef, (1, ∞))
        @test_throws ArgumentError Matrix{Float64}(undef, (∞, 1))

        @test Array{Float64}(undef, ()) isa Array{Float64,0}
        @test Array{Float64,0}(undef, ()) isa Array{Float64,0}
    end

    @testset "similar" begin
        a = 1:∞
        @test similar(a) isa CachedArray{Int}
        @test similar(a, Float64) isa CachedArray{Float64}
        @test similar(a, 5) isa Vector{Int}
        @test similar(a, (6,)) isa Vector{Int}
        @test similar(a, Float64, 5) isa Vector{Float64}
        @test similar(a, Float64, (6,)) isa Vector{Float64}
        @test similar(a, Float64, Base.OneTo(5)) isa Vector{Float64}
        @test similar(a, Float64, (Base.OneTo(5),)) isa Vector{Float64}
        @test similar(a, ∞) isa CachedArray{Int}
        @test similar(a, (∞,)) isa CachedArray{Int}
        @test similar(a, Float64, ∞) isa CachedArray{Float64}
        @test similar(a, Float64, (∞,)) isa CachedArray{Float64}
        @test similar(a, Float64, (∞,∞)) isa CachedArray{Float64}
        @test similar(a, Float64, oneto(∞)) isa CachedArray{Float64}
        @test similar(a, Float64, (oneto(∞),)) isa CachedArray{Float64}
        @test similar(a, Float64, (oneto(∞),oneto(∞))) isa CachedArray{Float64}

        @test similar([1,2,3],Float64,()) isa Array{Float64,0}

        @test similar(a, Float64, (2,∞)) isa CachedArray{Float64}
        @test similar(a, Float64, (∞,2)) isa CachedArray{Float64}

        @test similar(Array{Float64}, (∞,)) isa CachedArray{Float64}
        @test similar(Array{Float64}, (∞,∞)) isa CachedArray{Float64}
        @test similar(Array{Float64}, (2,∞)) isa CachedArray{Float64}
        @test similar(Array{Float64}, (∞,2)) isa CachedArray{Float64}

        @test similar(Array{Float64}, (oneto(∞),)) isa CachedArray{Float64}
        @test similar(Array{Float64}, (oneto(∞),oneto(∞))) isa CachedArray{Float64}
        @test similar(Array{Float64}, (oneto(∞),oneto(5))) isa CachedArray{Float64}
        @test similar(Array{Float64}, (oneto(5),oneto(∞))) isa CachedArray{Float64}
    end

    @testset "zeros/fill/ones" begin
        a = zeros(1,∞)
        @test length(a) === ℵ₀
        @test size(a) === (1,ℵ₀)
        @test a isa CachedArray{Float64}
        @test all(iszero,a[1,1:100])
        a[5] = 1
        @test a[1,1:100] == [zeros(4); 1; zeros(95)]

        a = fill(1,∞)
        @test length(a) === ℵ₀
        @test size(a) === (ℵ₀,)
        @test a isa CachedArray{Int}
        @test all(x -> x===1,a[1:100])
        a[5] = 2
        @test a[1:100] == [fill(1,4); 2; fill(1,95)]

        a = ones(∞)
        @test a isa CachedArray{Float64}
        a[5] = 2
        @test a[1:100] == [fill(1,4); 2; fill(1,95)]

        @test ones(5,∞)[:,1:10] == ones(5,10)
        @test ones(∞,5)[1:10,:] == ones(10,5)
        @test ones(∞,∞)[1:5,1:5] == ones(5,5)

        @test zeros(∞, 5)[1:10,:] == zeros(10,5)
        @test zeros(Int, ∞, 5)[1:10,:] == zeros(10,5)
        @test zeros(Int, 5, ∞)[:,1:10] == zeros(5, 10)
    end
end

@testset "ranges" begin
    @test size(10:1:∞) == (ℵ₀,)
    @testset "colon" begin
        @test @inferred(10:1:∞) === @inferred(range(10; step=1, length=ℵ₀))
        @inferred(1:0.2:∞)  === @inferred(range(1; step=0.2, length=ℵ₀))
        @inferred(1.0:0.2:∞)  === @inferred(range(1.0; step=0.2, length=ℵ₀))
        @inferred(1:∞) === @inferred(range(1; length=ℵ₀))
    end
    @test_throws ArgumentError 2:-.2:∞
    @test_throws ArgumentError 0.0:-∞
    @test_throws ArgumentError ∞:-1:1

    @test_throws ArgumentError (2:.2:-∞)
    @test_throws ArgumentError (2.:.2:-∞)
    @test_throws ArgumentError (1:-∞)

    @test ∞:1 ≡ 1:0

    @testset "indexing" begin
        @testset "axes" begin
            r = axes(big(1):∞,1)
            @test r == axes(r,1)
            @test r[typemax(Int)+big(1)] == typemax(Int)+big(1)
        end
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

        @test (1:∞)[Base.Slice(1:∞)] ≡ 1:∞
        @test Base.Slice(1:∞)[2:∞] ≡ 2:∞

        v = InfiniteArrays.InfStepRange(InfiniteCardinal{0}(), InfiniteCardinal{0}())
        @test v[1] == v[2] == InfiniteCardinal{0}()
    end
    @testset "length" begin
        @test length(.1:.1:∞) == ℵ₀
        @test length(1.1:1.1:∞) == ℵ₀
        @test length(1.1:1.3:∞) == ℵ₀
        @test length(1:1:∞) == ℵ₀
        @test length(1:.2:∞) == ℵ₀
        @test length(1.:.2:∞) == ℵ₀
        @test length(2:-.2:-∞) == ℵ₀
        @test length(2.:-.2:-∞) == ℵ₀

        @test Base.checked_length(1:∞) == length(1:∞)

        @testset "IteratorSize" begin
            @test Base.IteratorSize(1:2:∞) == Base.IsInfinite()
            @test Base.IteratorSize(1:∞) == Base.IsInfinite()
            s = Iterators.Stateful(2:∞)
            @test first(s) == 2
            @test first(s) == 3
        end
    end

    @testset "first" begin
        @test first(1:4:∞) == 1
        @test first(1:4:∞, 4) == range(1, step=4, length=4)
    end

    @testset "intersect" begin
        @test intersect(oneto(∞), 2:3) == intersect(2:3, oneto(∞)) == 2:3

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
        @test_throws ArgumentError intersect(1:6:∞, 0:4:∞) # supporting empty would break type inference

        @test intersect(1:∞,3) == 3:3
        @test intersect(1:∞, 2:∞, UnitRange(3,7), UnitRange(4,6)) == UnitRange(4,6)

        @test intersect(1:∞, 2) === intersect(2, 1:∞) === 2:2
        @test_skip intersect(1.0:∞, 2) == intersect(2, 1.0:∞) == [2.0] # gives infinite loop
    end

    @testset "sort/sort!/partialsort" begin
        @test sort(1:∞) ≡ sort!(1:∞) ≡ 1:∞
        @test sort(2:2:∞) ≡ sort!(2:2:∞) ≡ 2:2:∞
        @test_throws ArgumentError sort(2:-2:-∞)
        @test_throws ArgumentError sort!(2:-2:-∞)

        @testset "RangeCumsum" begin
            r = InfiniteArrays.OneToInf()
            rs = cumsum(r)
            @test sort(rs) === sort!(rs) === rs
            @test @inferred((rs -> Val(issorted(rs)))(rs)) isa Val{true}
            @test rs[end] ≡ ℵ₀
        end
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
        @test !(Complex(1, 1) in 1:∞)
        @test !(Complex(1, 1) in 1.0:∞)
        @test !(Complex(1.0, 1.0) in 1:∞)
        @test !(Complex(1.0, 1.0) in 1.0:∞)
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
        @test sum(1:∞) ≡ mean(1:∞) ≡ median(1:∞) ≡ ℵ₀
        @test sum(0:∞) ≡ mean(1:∞) ≡ median(1:∞) ≡ ℵ₀
        @test sum(0:2:∞) ≡ mean(0:2:∞) ≡ median(0:2:∞) ≡ RealInfinity(∞)
        @test sum(0:-2:-∞) ≡ mean(0:-2:-∞) ≡ median(0:-2:-∞) ≡ -∞
    end

    @testset "broadcasted operations with scalars" begin
        @test Base.BroadcastStyle(typeof(1:∞)) isa LazyArrayStyle{1}
        @test Base.BroadcastStyle(typeof(Base.Slice(1:∞))) isa LazyArrayStyle{1}

        @test broadcast(-, 1:∞, 2) ≡ -1:∞
        @test broadcast(-, 1:∞, 0.25) ≡ 1-0.25:∞
        @test broadcast(+, 1:∞, 2) ≡ 3:∞
        @test broadcast(+, 1:∞, 0.25) ≡ 1+0.25:∞
        @test broadcast(+, 1:2:∞, 1) ≡ 2:2:∞
        @test broadcast(+, 1:2:∞, 0.3) ≡ 1+0.3:2:∞
        @test broadcast(-, 1:2:∞, 1) ≡ 0:2:∞
        @test broadcast(-, 1:2:∞, 0.3) ≡ 1-0.3:2:∞
        @test broadcast(-, 2, 1:∞) ≡ 1:-1:-∞
        @test exp.((1:∞)') ≡ broadcast(exp, (1:∞)') ≡ exp.(1:∞)'
        @test exp.(transpose(1:∞)) ≡ broadcast(exp, transpose(1:∞)) ≡ transpose(exp.(1:∞))
        @test 1 .+ (1:∞)' ≡ broadcast(+, 1, (1:∞)') ≡ (2:∞)'
        @test 1 .+ transpose(1:∞) ≡ broadcast(+, 1, transpose(1:∞)) ≡ transpose(2:∞)
        @test (1:∞)' .+ 1 ≡ broadcast(+, (1:∞)', 1) ≡ (2:∞)'
        @test transpose(1:∞) .+ 1 ≡ broadcast(+, transpose(1:∞), 1) ≡ transpose(2:∞)
    end

    @testset "near-equal ranges" begin
        @test 0.0:0.1:∞ != 0.0f0:0.1f0:∞
    end

    @testset "constprop in comparing OneToInf" begin
        r1 = OneToInf{Int8}()
        r2 = OneToInf{Int16}()
        v = @inferred ((r1,r2) -> Val(r1 == r2))(r1, r2)
        @test v isa Val{true}
    end

    @testset "comparing InfiniteUnitRanges and OneToInf" begin
        @test 1:2:∞ == 1:2:∞ != 1:3:∞ != 2:3:∞ == 2:3:∞ != 2:∞
        @test 1:1:∞ == 1:∞ == 1:∞ == OneToInf() == OneToInf()
    end

    @testset "Base.OneTo (misleading) overrides" begin
        @test_skip Base.OneTo{BigInt}(∞) isa OneToInf{BigInt}
        @test oneto(∞) isa OneToInf{Int}
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

        @test AbstractArray{Float64}(1:2:∞) ≡ AbstractVector{Float64}(1:2:∞) ≡
                convert(AbstractVector{Float64}, 1:2:∞) ≡ convert(AbstractArray{Float64}, 1:2:∞)

        @test unitrange(oneto(∞)) ≡ InfUnitRange(oneto(∞)) ≡ InfUnitRange{Int}(oneto(∞)) ≡ InfUnitRange(1)

        @test oneto(RealInfinity()) ≡ oneto(ComplexInfinity()) ≡ OneToInf()
        @test_throws ArgumentError oneto(-ComplexInfinity())
        @test_throws ArgumentError oneto(-∞)
    end

    @testset "inf-range[inf-range]" begin
        @test (1:∞)[1:∞] == 1:∞
        @test (1:∞)[OneToInf()] == 1:∞
        @test (1:∞)[2:∞] == 2:∞
        @test_throws BoundsError (1:∞)[-1:∞]
        @test (1:-1:-∞)[1:∞] == 1:-1:-∞
    end

    @testset "view(InfStepRange, inf-range)" begin
        r = 2:5:∞
        @test view(r, axes(r,1)) === r
        @test view(r, 1:1:∞) === r
        @test view(r, 2:2:∞) === 7:10:∞
        r2 = 2:5:10_000 # arbitrary high upper cutoff
        @test view(r, 4:10) == view(r2, 4:10)
        @test view(r, 4:7:50) == view(r2, 4:7:50)
    end

    @testset "OneToInf" begin
        r = OneToInf()
        @test !isempty(r)
        @test length(r) == ℵ₀
        @test size(r) == (ℵ₀,)
        @test step(r) == 1
        @test first(r) == 1
        @test last(r) == ∞
        @test minimum(r) == 1
        @test maximum(r) == ℵ₀
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

        @test Base.unsafe_indices(Base.Slice(r)) == (r,)

        @testset "iteration with zip + finite iterator" begin
            z = zip(OneToInf(), 1:100)
            @test axes(z) == (Base.OneTo(100),)
            @test size(z) == (100,)
        end
    end

    @testset "show" begin
        # NOTE: Interpolating Int to ensure it's displayed properly across 32- and 64-bit
        @test summary(1:∞) == "ℵ₀-element InfUnitRange{$Int} with indices OneToInf()"
        @test Base.inds2string(axes(1:∞)) == "OneToInf()"
    end

    @testset "end" begin
        @test oneto(∞)[end] ≡ oneto(∞)[∞] ≡ ℵ₀
        @test (1:∞)[end] ≡ (1:∞)[∞] ≡ ℵ₀
        @test (1:2:∞)[end] ≡ (1:2:∞)[∞] ≡ ℵ₀
        @test (1.0:2:∞)[end] ≡ (1.0:2:∞)[∞] ≡ ℵ₀
    end

    @testset "union" begin
        @test @inferred((1:∞) ∪ (3:∞)) ≡ @inferred((3:∞) ∪ (1:∞)) ≡ 1:∞
        @test @inferred((1:∞) ∪ (3:1:∞)) ≡ @inferred((3:1:∞) ∪ (1:∞)) ≡ 1:1:∞
        @test @inferred((2:2:∞) ∪ (4:2:∞)) ≡ 2:2:∞
        @test (2.0:1.5:∞) ∪ (3.5:1.5:∞) ≡ 2.0:1.5:∞
        @test (1:∞) ∪ (2:2:∞) ≡ 1:1:∞
        @test (6:4:∞) ∪ (2:2:∞) ≡ 2:2:∞
        @test_throws ArgumentError (3:∞) ∪ (2:2:∞)
        @test_throws ArgumentError (2:2:∞) ∪ (3:∞)
        @test_throws ArgumentError (2:3:∞) ∪ (2:2:∞)

        @test @inferred(union(1:∞)) ≡ 1:∞
        @test @inferred(union(1.5:∞)) ≡ 1.5:∞
        @test @inferred(union(1:∞, 2:∞, 4:∞)) ≡ 1:∞
        @test @inferred(union(1.5:∞, 0.5:∞, 2.5:2:∞)) ≡ 0.5:1:∞
    end

    @testset "adjoint indexing" begin
        a = (1:∞)'
        @test a[:,:] ≡ a
        @test a[1,:] ≡ 1:∞
        @test a[1,2:2:end] ≡ 2:2:∞
        @test a[:,2:∞][1,1:10] == a[2:11]
    end

    @testset "big" begin
        @test (big(1):∞)[5] isa BigInt
        @test range(big(1), ∞; step=1) == big(1):1:∞
        @test range(big(1.0), ∞; step=2.5) == range(big(1.0); step=big(2.5), length=ℵ₀) == range(big(1.0), ∞; step=big(2.5)) == big(1.0):2.5:∞
    end

    @testset "maximum/minimum" begin
        @test maximum(2:∞) ≡ ℵ₀
        @test minimum(2:∞) ≡ 2
    end

    @testset "getindex[∞]" begin
        @test_throws BoundsError (2:6)[∞]
        @test (2:∞)[∞] ≡ ℵ₀
        @test oneto(∞)[∞] ≡ ℵ₀
        @test oneto(∞)[RealInfinity()] ≡ ℵ₀
    end

    @testset "show" begin
        @test stringmime("text/plain", 2:∞) ≡ "2:∞"
        @test stringmime("text/plain", OneToInf{BigInt}()) ≡ "OneToInf{BigInt}()"
    end

    @testset "in" begin
        @test ∞ ∉ (1:∞)
    end

    @testset "iterate" begin
        for r in (oneto(∞), 1:∞, 1:1:∞)
            x = 0
            for k in r
                x += 1
                k > 5 && break
            end
            @test x == 6
        end
    end

    @testset "vcat" begin
        @test [1:∞;] === 1:∞

        @testset for r in (2, 2.0)
            v = [r; 1:∞]
            @test v isa AbstractVector{typeof(r)}
            @test isinf(length(v))
            @test v[1] == r
            if typeof(v) == Int # fast infinite getindex is not defined for Float64
                @test v[2:∞] == 1:∞
            else
                @test v[2:10] == 1:9
            end
        end

        @testset for r in (1:2, [1,2])
            v = [r; 1:∞]
            @test v isa AbstractVector{Int}
            @test isinf(length(v))
            @test v[axes(r,1)] == r
            if isfinite(length(r))
                @test v[length(r) .+ 1:∞] == 1:∞
            end
        end

        @testset for r in (1.0:2.0, [1.0,2.0])
            v = [r; 1:∞]
            @test v isa AbstractVector{Float64}
            @test isinf(length(v))
            if isfinite(length(r))
                @test v[axes(r,1)] == r
            end
            @test v[length(r) .+ (1:10)] == 1:10
        end

        @test_throws ArgumentError [1:∞; 1]
        @test_throws ArgumentError [1:∞; 1:∞]
        @test_throws ArgumentError [1:∞; 1]
        @test_throws ArgumentError [1:∞; 1:∞; 1:∞]
        @test_throws ArgumentError [1:∞; 1:∞; 1:∞; 1]
        @test_throws ArgumentError [1:∞; 1.0:∞]
        @test_throws ArgumentError [1:∞; 1:2]
        @test_throws ArgumentError [1:∞; 1:2; 1]
        @test_throws ArgumentError [1:∞; 1:2; 1:∞]
        @test_throws ArgumentError [1:∞; 1:2.0]
        @test_throws ArgumentError [1:∞; 1:2.0; 1:∞]
        @test_throws ArgumentError [1:∞; 1:2.0; 1]
        @test_throws ArgumentError [1:∞; [1]]
        @test_throws ArgumentError [1:∞; [1]; 1:∞]
        @test_throws ArgumentError [1:∞; [1.0]]
        @test_throws ArgumentError [1:∞; [1.0]; 1:∞]
        @test_throws ArgumentError [1:∞; 1; [1]]
    end
end

@testset "fill" begin
    @testset "fill sizes" begin
        for A in (Zeros(∞), Fill(1,∞), Ones(∞),
                    Zeros(5,∞), Ones(5,∞), Fill(1,5,∞),
                    Zeros(∞,5), Ones(∞,5), Fill(1,∞,5),
                    Zeros(∞,∞), Ones(∞,∞), Fill(1,∞,∞))
            @test length(A) ≡ ℵ₀
        end
        @test size(Zeros(∞,5)) ≡ (ℵ₀,5)
        @test size(Zeros(5,∞)) ≡ (5,ℵ₀)
    end

    @testset "Fill indexing" begin
        B = Ones(∞,∞)
        @test IndexStyle(B) == IndexCartesian()
        V = view(B,:,1)
        @test_broken size(V) == (ℵ₀,1)
        V = view(B,1,:)
        @test size(V) == (ℵ₀,)
        V = view(B,1:1,:)
        @test size(V) == (1,ℵ₀)
    end

    @testset "Fill reindex" begin
        F = Fill(2.0,2,∞)
        @test reshape(F,∞) ≡ reshape(F,OneToInf()) ≡ reshape(F,(OneToInf(),)) ≡ reshape(F,Val(1)) ≡ Fill(2.0,∞)
    end

    @testset "adjtrans copy" begin
        @test copy((1:∞)') ≡ (1:∞)'
        @test copy(transpose(1:∞)) ≡ transpose(1:∞)
    end

    @testset "Fill slices" begin
        A = Fill(2,∞,∞)
        Z = Zeros(∞,∞)
        @test A[:,1] ≡ A[1,:] ≡ A[1:∞,1] ≡ Fill(2,∞)
        @test Z[:,1] ≡ Z[1,:] ≡ Z[1:∞,1] ≡ Zeros(∞)

        @test A[2:∞,1:∞] ≡ A[2:∞,:] ≡ A[:,1:∞] ≡  A
        @test A[5,2:∞] ≡ A[2:∞,5] ≡ Fill(2,∞)
    end

    @testset "maximum/minimum/sum" begin
        c = cache(Fill(2,∞));
        c[1] = 1;
        @test maximum(c) == maximum(Vcat([1], Fill(2,∞))) == 2
        c[1:3] = 1:3;
        @test maximum(c) == maximum(Vcat([1,2,3], Fill(2,∞))) == 3
        @test minimum(c) == minimum(Vcat([1,2,3], Fill(2,∞))) == 1
        @test sum([1; zeros(∞)]) ≡ 1.0
        @test sum([1; ones(∞)]) ≡ 1.0∞
    end
end

@testset "diagonal" begin
    D = Diagonal(1:∞)
    @test D[1:10,1:10] == Diagonal(1:10)
    @test D[:,1:5][2:5,:] == D[2:5,1:5]
    @test D[1:5,:][:,2:5] == D[1:5,2:5]
    @test D[:,:][1:5,1:5] == D[1:5,1:5]

    @test D[:,5][1:10] == D[1:10,5]
    @test D[5,:][1:10] == D[5,1:10]

    @test D^2 isa Diagonal
    @test D*D isa Diagonal
    @test MemoryLayout(typeof(D.diag)) == LazyLayout()
    @test MemoryLayout(typeof(D)) == DiagonalLayout{LazyLayout}()
    @test Base.BroadcastStyle(typeof(D)) == LazyArrayStyle{2}()
    @test Base.BroadcastStyle(typeof(permutedims(D.diag))) == LazyArrayStyle{2}()
    bc = broadcasted(*,Ones(∞,∞),permutedims(D.diag))
    @test bc isa Broadcasted{LazyArrayStyle{2}}
    @test instantiate(bc) isa Broadcasted{LazyArrayStyle{2}}
    @test copy(instantiate(bc)) isa BroadcastArray
    @test broadcast(*,Ones(∞,∞),permutedims(D.diag)) isa BroadcastArray
    @test Ones(∞,∞)*D isa BroadcastArray
    @test (Ones(∞,∞)*D)[1:10,1:10] == Ones(10,10)*D[1:10,1:10]
    @test @inferred(broadcast(*,Ones{Int}(∞),D)) ≡ @inferred(broadcast(*,D,Ones{Int}(∞))) ≡ D
    @test @inferred(broadcast(*,Ones(∞),D)) == @inferred(broadcast(*,D,Ones(∞))) == Diagonal(1.0:∞)
    @test @inferred(broadcast(*,Ones{Int}(∞)',D)) == @inferred(broadcast(*,D,Ones{Int}(∞)')) == D
    @test @inferred(broadcast(*,Ones(∞)',D)) == @inferred(broadcast(*,D,Ones(∞)')) == Diagonal(1.0:∞)
    @test @inferred(broadcast(*,Fill(2,∞)',D)) ≡ @inferred(broadcast(*,D,Fill(2,∞)')) ≡ 2D

    @test Eye{Int}(∞, ∞) isa Eye{Int}
    @test Eye{Int}(∞, 5) isa Eye{Int}
    @test Eye{Int}(5, ∞) isa Eye{Int}
    @test Eye(∞, ∞) isa Eye{Float64}
    @test Eye(∞, 5) isa Eye{Float64}
    @test Eye(5, ∞) isa Eye{Float64}

    @test Eye{Int}(∞) * D ≡ Eye{Int}(∞) * D ≡ D
    @test Eye(∞) * D == Eye(∞) * D == D
    @test Eye(∞) == Eye(∞)^0 == Eye(∞)^1 == Eye(∞)^2 == one(Eye(∞)) == copy(Eye(∞)) == one(Diagonal(Fill(2,∞)))
    @test Diagonal(Fill(2,∞)) == copy(Diagonal(Fill(2,∞)))


    @test permutedims(D) ≡ D
    @test copy(D) ≡ D

    @test 2D ≡ D*2 ≡ 2 .* D ≡ D .* 2
end

@testset "concat" begin
    @testset "concat indexing" begin
        A = Vcat(1:10, 1:∞)
        @test @inferred(length(A)) == ℵ₀
        @test @inferred(A[5]) == A[15] == 5
        @test A[end] == @inferred(A[∞]) == ∞
        @test_throws BoundsError Vcat(1:10)[∞]

        A = Vcat(Ones(1,∞), Zeros(2,∞))
        @test @inferred(size(A)) == (3,ℵ₀)
        @test @inferred(A[1,5]) == 1
        @test @inferred(A[3,5]) == 0
        @test_throws BoundsError A[4,1]

        A = Vcat(Ones{Int}(1,∞), Diagonal(1:∞))
        @test @inferred(size(A)) ≡ (ℵ₀,ℵ₀)
        @test @inferred(A[1,5]) ≡ 1
        @test @inferred(A[5,5]) ≡ 0
        @test @inferred(A[6,5]) ≡ 5
        @test_throws BoundsError A[-1,1]

        A = Vcat(Ones{Float64}(1,∞), Diagonal(1:∞))
        @test @inferred(size(A)) ≡ (ℵ₀,ℵ₀)
        @test @inferred(A[1,5]) ≡ 1.0
        @test @inferred(A[5,5]) ≡ 0.0
        @test @inferred(A[6,5]) ≡ 5.0
        @test_throws BoundsError A[-1,1]

        A = Vcat(1, Zeros(∞))
        @test @inferred(A[1]) ≡ 1.0
        @test @inferred(A[2]) ≡ 0.0

        A = Hcat(Ones(∞), Zeros(∞,2))
        @test @inferred(size(A)) == (ℵ₀,3)
        @test @inferred(A[5,1]) == 1
        @test @inferred(A[5,3]) == 0
        @test_throws BoundsError A[1,4]

        A = Hcat(Ones{Int}(∞), Diagonal(1:∞))
        @test @inferred(size(A)) ≡ (ℵ₀,ℵ₀)
        @test @inferred(A[5,1]) ≡ 1
        @test @inferred(A[5,5]) ≡ 0
        @test @inferred(A[5,6]) ≡ 5
        @test_throws BoundsError A[-1,1]

        A = Hcat(1, (1:∞)')
        @test A[1,:] isa Vcat{<:Any,1}
        @test A[1,:][1:10] == A[1,1:10]
    end

    # This should be generalized, but it at the moment
    # it is restricted to a single Number. Support smart
    # addition for any number of Number/SVector's would be better
    # allowibng for the tail to be variable length
    @testset "Vcat special case" begin
        @test Vcat(1,Zeros{Int}(∞)) + Vcat(3,Zeros{Int}(∞)) ≡
            Vcat(1,Zeros{Int}(∞)) .+ Vcat(3,Zeros{Int}(∞)) ≡
            Vcat(4,Zeros{Int}(∞))

        @test size(Vcat(1:∞)) ≡ (ℵ₀,)
    end

    @testset "Vcat infrange getindex" begin
        x = Vcat(1, Fill(2,∞))
        @test x[1:end] ≡ x[1:∞] ≡ x
        @test x[3:end] ≡ x[3:∞] ≡ Fill(2,∞)
    end

    @testset "maximum/minimum Vcat" begin
        x = Vcat(1:2, [1,1,1,1,1], 3, Fill(4,∞))
        @test maximum(x) == 4
        @test minimum(x) == 1
    end

    @testset "special vcat" begin
        @test [1; Zeros(∞)][1:10] == [1; zeros(9)]
        @test [[1,2,3]; Zeros(∞)][1:10] == [1;2;3;zeros(7)]
        @test [1; zeros(∞)] isa CachedArray
        @test [[1,2,3]; zeros(∞)] isa CachedArray

        @test [1; 2; zeros(Int,∞)] isa CachedArray
        @test [1; 2; 3; zeros(Int,∞)] isa CachedArray
        @test [[1,2]; 3; zeros(Int,∞)] isa CachedArray
        @test [2; [1,2]; 3; zeros(Int,∞)] isa CachedArray

        @test [[1,2]; [3,4]; Zeros(∞)] isa Vcat{<:Any,1,<:Tuple{Array,Zeros}}
        @test [[1,2]; [3,4]; [5,6]; Zeros(∞)] isa Vcat{<:Any,1,<:Tuple{Array,Zeros}}

        @test [randn(2,2); Zeros(∞,2)] isa Vcat{<:Any,2,<:Tuple{Array,Zeros}}

        a = [[1,2,3]; zeros(Int,∞)]
        @test a[3:∞][1:5] == a[3:7]
        @test cache(1:∞)[2:∞][1:5] == 2:6

        D = Diagonal(1:∞)
        @test [D[2:5,:]; D][1:10,1:10] == [D[2:5,1:10]; D[1:6,1:10]]
        @test [D[3,:] D][1:10,1:10] == [D[3,1:10] D[1:10,1:9]]
        @test [D[:,1:2] D][1:10,1:10] == [D[1:10,1:2] D[1:10,1:8]]

        @test [1:5 D[1:5,:]][:,1:5] == [1:5 D[1:5,1:4]]

        @test [1:∞ D[:,1:5]][1:10,:] == [1:10 D[1:10,1:5]]
        @test [1 Zeros(1,∞)][:,1:10] == [1 zeros(1,9)]

        @test [[1; zeros(∞)] D[:,1:5]][1:10,:] == [[1; zeros(9)] D[1:10,1:5]]
        @test [[1; zeros(∞)] BandedMatrix(D[:,1:5])][1:10,:] == [[1; zeros(9)] D[1:10,1:5]]

        @test cat([1,2,3],zeros(∞); dims=1) == cat(1:3,zeros(∞); dims=1) == [[1,2,3]; zeros(∞)]
    end

    @testset "sparse print" begin
        A = Vcat(1, Zeros(∞))
        @test colsupport(A,1) == 1:1
        @test Base.replace_in_print_matrix(A, 2, 1, "0") == "⋅"
        @test stringmime("text/plain", A; context=(:limit => true)) ==
                "vcat($Int, ℵ₀-element Zeros{Float64, 1, Tuple{OneToInf{$Int}}} with indices OneToInf()) with indices OneToInf():\n 1.0\n  ⋅ \n  ⋅ \n  ⋅ \n  ⋅ \n  ⋅ \n  ⋅ \n  ⋅ \n  ⋅ \n  ⋅ \n ⋮"
        A = Vcat(Ones{Int}(1,∞), Diagonal(1:∞))
        @test Base.replace_in_print_matrix(A, 2, 2, "0") == "⋅"
    end

    @testset "copymutable" begin
        @test Base.copymutable(Vcat(1., Zeros(∞))) isa CachedArray
        @test Base.copymutable(Vcat([1.], Zeros(∞))) isa CachedArray
        @test Base.copymutable(Vcat([1.,2.], zeros(∞))) isa CachedArray
        @test Base.copymutable(Vcat(1.,2., zeros(∞))) isa CachedArray
    end

    @testset "infinite indexing" begin
        a = Vcat(1, 1:∞)
        @test a[:] isa Vcat
        @test a[3:∞] ≡ 2:∞
        @test a[3:2:∞] isa Vcat

        A = Vcat(Ones(1,∞), Fill(2,1,∞))
        @test A[:,:] == A
        @test A[:,2:∞] isa Vcat

        A = Vcat(Ones(5,5), Fill(2,∞,5))
        @test A[:,:] == A
        @test A[2:∞,:] isa Vcat

        A = Vcat(Ones(1,∞), Fill(2,∞,∞))
        @test A[:,:] == A
        @test A[2:∞,2:∞] isa Vcat
        @test A[2:∞,2:∞][1:10,1:10] == fill(2,10,10)

        B = Hcat(Ones(∞), Diagonal(1:∞))
        @test B[2:∞,1:∞][1:10,1:10] == B[2:11,1:10]
        @test B[2:5,1:∞][:,1:10] == B[2:5,1:10]
        @test B[2:∞,2:5][1:10,:] == B[2:11,2:5]

        @test B[3,1:∞][1:10] == B[3,1:10]
        @test B[1:∞,3][1:10] == B[1:10,3]
    end

    @testset "adjoint copy"  begin
        a = Vcat(1,(1:∞))'
        b = transpose(Vcat(1,(1:∞)))
        @test copy(a) ≡ a
        @test copy(b) ≡ b
    end

    @testset "hcat" begin
        @test [Zeros(∞) Diagonal(1:∞)][1:10,1:11] == [zeros(10) Diagonal(1:10)]
        @test [Zeros(∞,3) Diagonal(1:∞)][1:10,1:13] == [zeros(10,3) Diagonal(1:10)]
        @test [Fill(2,∞) Diagonal(1:∞)][1:10,1:11] == [fill(2,10) Diagonal(1:10)]
        @test [Fill(2,∞,3) Diagonal(1:∞)][1:10,1:13] == [fill(2,10,3) Diagonal(1:10)]
    end

    @testset "Banded concat" begin
        A = _BandedMatrix((0:∞)', ℵ₀, -1, 1)
        D = Diagonal(1:∞)
        @test [view(D,1:1,:); A][1:10,1:10] == [D[1:1,1:10]; A[1:9,1:10]]
        @test [view(A,1:1,:); A][1:10,1:10] == [A[1:1,1:10]; A[1:9,1:10]]
        @test [Ones(∞) A][1:10,1:10] == [ones(10) A[1:10,1:9]]
        @test [Ones(∞,2) A][1:10,1:10] == [ones(10,2) A[1:10,1:8]]
        @test [view(D,1,:) A][1:10,1:10] == [D[1,1:10] A[1:10,1:9]]
        @test [view(A,1,:) A][1:10,1:10] == [A[1,1:10] A[1:10,1:9]]
        @test [view(D,:,1:1) A][1:10,1:10] == [D[1:10,1] A[1:10,1:9]]
        @test [view(A,:,1:1) A][1:10,1:10] == [A[1:10,1] A[1:10,1:9]]
        @test [view(A,1,:) A][1:10,1:10] == [A[1,1:10] A[1:10,1:9]]
    end
end

@testset "broadcasting" begin
    @testset "∞ BroadcastArray" begin
        A = 1:∞
        B = BroadcastArray(exp, A)
        @test length(B) == ℵ₀
        @test B[6] == exp(6)
        @test exp.(A) ≡ B
        @test B[2:∞] isa BroadcastArray
        B = Diagonal(1:∞) .+ 1
        @test B isa BroadcastArray{Int}
        @test B[1,5] ≡ 1
        @test B[6,6] == 6+1
        B = Diagonal(1:∞) - Ones{Int}(∞,∞) # lowers to broadcast
        @test B isa BroadcastArray{Int}
        @test B[1,5] ≡ -1
        @test B[6,6] == 6-1
    end

    @testset "Broadcast Fill Lowers" begin
        @test broadcast(+, Zeros{Int}(∞) , Fill(1,∞)) isa Fill
        @test broadcast(+, Zeros{Int}(∞) , Zeros(∞)) isa Zeros
        @test broadcast(*, Ones(∞), Ones(∞)) ≡ Ones(∞)
        @test broadcast(*, Ones{Int}(∞), 1:∞) ≡ broadcast(*, 1:∞, Ones{Int}(∞)) ≡ 1:∞
        @test broadcast(*, Fill(2,∞), 1:∞) ≡ broadcast(*, 1:∞, Fill(2,∞)) ≡ 2:2:∞
        @test broadcast(*, Fill([1,2],∞), 1:∞) isa BroadcastVector
        @test broadcast(*, Fill([1,2],∞), 1:∞)[1:3] == broadcast(*, 1:∞, Fill([1,2],∞))[1:3] == [[1,2],[2,4],[3,6]]

        @test broadcast(*, 1:∞, Ones(∞)') isa BroadcastArray
        @test broadcast(*, 1:∞, Fill(2,∞)') isa BroadcastArray
        @test broadcast(*, Diagonal(1:∞), Ones{Int}(∞)') ≡ broadcast(*, Ones{Int}(∞)', Diagonal(1:∞)) ≡ Diagonal(1:∞)
        @test broadcast(*, Diagonal(1:∞), Fill(2,∞)') ≡ broadcast(*, Fill(2,∞)', Diagonal(1:∞)) ≡ Diagonal(2:2:∞)
    end

    @testset "subview inf broadcast" begin
        b = BroadcastArray(exp, 1:∞)
        v = view(b, 3:∞) .+ 1
        @test v isa BroadcastArray
        @test b[3:10] .+ 1 == v[1:8]
    end

    @testset "views of matrices" begin
        D = Diagonal(1:∞)
        V = Vcat(Ones(2,∞), D)
        @test view(D,:,5) .+ 1 isa BroadcastVector
        @test view(D,5,:) .+ 1 isa BroadcastVector
        @test view(V,:,5) .+ 1 isa BroadcastVector
        @test view(V,5,:) .+ 1 isa BroadcastVector

        @test view(D,2:∞,2:∞) .+ 1 isa BroadcastMatrix
        @test view(V,2:∞,2:∞) .+ 1 isa BroadcastMatrix

        @test view(D,2:∞,[1,2,3]) .+ 1 isa BroadcastMatrix
        @test view(D,[1,2,3],2:∞) .+ 1 isa BroadcastMatrix
        @test view(V,2:∞,[1,2,3]) .+ 1 isa BroadcastMatrix
        @test view(V,[1,2,3],2:∞) .+ 1 isa BroadcastMatrix
    end

    @testset "inf broadcast views" begin
        a = BroadcastArray(cos, 1:∞)
        r = div.(1:∞, 2) .+ 1
        b = SubArray(a, (r,))
        @test b[1:6] == a[r[1:6]]
        @test_broken Base.BroadcastStyle(typeof(b)) isa LazyArrayStyle
        c = SubArray(a, (view(r,2:∞),))
        @test c[1:6] == a[r[2:7]]
        @test Base.BroadcastStyle(typeof(c)) isa LazyArrayStyle
    end

    @testset "structured matrices" begin
        r = 1:∞
        f = Fill(2, ∞)
        for B in (Bidiagonal(r, r, :U), Tridiagonal(r, r, r), SymTridiagonal(r, r),
                   Bidiagonal(f, f, :U), Tridiagonal(f, f, f), SymTridiagonal(f, f))
            B2 = B .+ B
            @test B2[1:10, 1:10] == 2B[1:10, 1:10]
        end
    end
end

@testset "Cumsum and diff" begin
    @test cumsum(Ones(∞)) ≡ 1.0:1.0:∞
    @test cumsum(Fill(2,∞)) ≡ 2:2:∞
    @test cumsum(Ones{Int}(∞)) ≡ oneto(∞)
    @test cumsum(Ones{BigInt}(∞)) ≡ OneToInf{BigInt}()

    @test diff(oneto(∞)) ≡ Ones{Int}(∞)
    @test diff(1:∞) ≡ Ones{Int}(∞)
    @test diff(1:2:∞) ≡ Fill(2,∞)
    @test diff(1:2.0:∞) ≡ Fill(2.0,∞)

    x = Vcat([3,4], Ones{Int}(5), 3, Fill(2,∞))
    y = @inferred(cumsum(x))
    @test y isa Vcat
    @test y[1:12] == cumsum(x[1:12])

    @test diff(x[1:10]) == diff(x)[1:9]
    @test diff(y)[1:20] == x[2:21]

    @test cumsum(x).args[2] ≡ 8:12
    @test last(y.args) == sum(x[1:9]):2:∞

    for r in (3:4:∞, 2:∞, oneto(∞))
        c = cumsum(r)
        @test c isa InfiniteArrays.RangeCumsum
        @test c[Base.OneTo(20)] == c[1:20] == [c[k] for k=1:20] == cumsum(r[1:20])
        @test c[2:20] == [c[k] for k=2:20] == cumsum(r[1:20])[2:end]
        @test c == c
        @test c[Base.OneTo(20)] isa InfiniteArrays.RangeCumsum
        @test exp.(c)[1:20] == exp.(c[1:20])
    end

    @test cumsum(3:4:∞)[end] ≡ cumsum(3:4:∞)[∞] ≡ cumsum(3:4:∞)[ℵ₀] ≡ RealInfinity()
    @test cumsum(2:∞)[end] ≡ cumsum(2:∞)[∞] ≡ cumsum(2:∞)[ℵ₀] ≡ cumsum(oneto(∞))[end] ≡ cumsum(oneto(∞))[∞] ≡ cumsum(oneto(∞))[ℵ₀] ≡  ℵ₀

    @test_throws BoundsError cumsum(oneto(∞))[-5]

    @test cumsum(1:∞)[2:∞][1:5] == cumsum(1:6)[2:end]

    @testset "union of cumsum" begin
        r1 = InfiniteArrays.OneToInf{Int8}()
        r2 = InfiniteArrays.OneToInf{Int16}()
        rs = union(cumsum(r1), cumsum(r2))
        @test rs == cumsum(r2)
    end
end

@testset "Sub-array" begin
    @test Ones(∞)[3:∞] ≡ Ones(∞)
    @test Ones{Int}(∞)[4:6] ≡ Ones{Int}(3)
    @test (1:∞)[3:∞] ≡ 3:∞
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

@testset "show" begin
    @test repr(Vcat(1:∞)) == "[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, …]"
    @test repr(Vcat(2,1:∞)) == "[2, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, …]"
end

@testset "*" begin
    A = Fill(1,3,∞)
    B = Diagonal(1:∞)
    C = Vcat([1,2,3], Zeros(∞))
    D = Vcat(Fill(1,3,∞), Zeros(∞,∞))

    AB = A*B
    @test AB isa BroadcastArray
    @test size(AB) == (3,ℵ₀)
    @test (AB)[1:3,1:10] == Fill(1,3,10)*Diagonal(1:10)

    @test size(A*B*C) == (3,)
    @test (A*B*C)[1] == 14
    @test A*B*C == fill(14,3)
    @test_throws BoundsError (A*B*C)[1:10,1:10]
    @test A*B*D isa ApplyArray
    @test (A*B*D)[1:3,1:5] == fill(6.0,3,5)
end

@testset "MemoryLayout" begin
    @test MemoryLayout(OneToInf{Int}) == LazyLayout()
    @test MemoryLayout(0:∞) == LazyLayout()
    @test MemoryLayout((0:∞)') == DualLayout{LazyLayout}()
    A = _BandedMatrix((0:∞)', ℵ₀, -1, 1)
    @test MemoryLayout(A) == BandedColumns{LazyLayout}()

    @test A[2:∞,1:∞][1:10,1:10] == A[2:11,1:10]

    @test A[2,1:∞][1:10] == A[2,1:10]
    @test A[1:∞,2][1:10] == A[1:10,2]

    @test A[[2,3],1:∞][:,1:10] == A[[2,3],1:10]
    @test A[1:∞,[2,3]][1:10,:] == A[1:10,[2,3]]
end

@testset "Banded" begin
    A = _BandedMatrix((0:∞)', ℵ₀, -1, 1)
    @test (Eye{Int}(∞) * A).data ≡ A.data
    @test 2.0A isa BandedMatrix
    @test (2.0A)[1:10,1:10] == 2.0A[1:10,1:10]
    @test 2.0\A isa BandedMatrix
    @test (2.0\A)[1:10,1:10] == 2.0\A[1:10,1:10]
    @test A/2 isa BandedMatrix
    @test (A/2)[1:10,1:10] == (A/2)[1:10,1:10]

    @test A * Eye(∞) isa BandedMatrix
    @test A / Eye(∞) isa BandedMatrix
    @test Eye(∞) * A isa BandedMatrix
    @test Eye(∞) \ A isa BandedMatrix
end

@testset "reshaped" begin
    @test InfiniteArrays.ReshapedArray(1:6,(2,3)) == [1 3 5; 2 4 6]
    @test InfiniteArrays.ReshapedArray(1:∞,(1,ℵ₀))[1,1:10] == 1:10
    @test reshape(1:∞,1,∞) ≡ InfiniteArrays.ReshapedArray(1:∞,(1,ℵ₀))
    @test parentindices(reshape(1:∞,1,∞)) ≡ (oneto(∞),)
    @test permutedims(1:∞) isa InfiniteArrays.ReshapedArray
    @test permutedims(1:∞)[1,1:10] == (1:10)
    a = reshape(Vcat(Fill(1,1,∞),Fill(2,2,∞)),∞)
    @test a[1:7] == [1, 2, 2, 1, 2, 2, 1]
    @test permutedims(permutedims(1:∞)) ≡ 1:∞
    @test parentindices(a) ≡ (oneto(3),oneto(∞))
    @test Base.unaliascopy(a) ≡ a
    @test Base.dataids(a) == Base.dataids(parent(a))
    @test a[Base.ReshapedIndex(5)] == a[5]

    b = reshape([1; zeros(∞)],1,∞)
    b[1,5] = 6
    @test parent(b)[5] == 6

    @test reshape(Ones(∞), 1, ∞) ≡ Ones(1,∞)
    @test reshape(Ones(∞), (1, ∞)) ≡ Ones(1,∞)
end

@testset "norm/dot" begin
    for p in (-Inf, 0, 0.1, 1, 2, 3, Inf)
        @test norm(Zeros(∞), p) == 0.0
        @test norm(Fill(5),p) ≈ norm(Array(Fill(5)),p) # tests tuple bug
        @test norm(Zeros{Float64}(),p) == 0.0 # tests tuple bug
    end
    @test norm([1; zeros(∞)]) ≡ 1.0
    @test dot([1; zeros(∞)], [1; zeros(∞)]) ≡ 1.0
    @test dot([1; zeros(∞)], 1:∞) ≡ 1.0
    @test dot(1:∞, [1; zeros(∞)]) ≡ 1.0
end

@testset "sub-Eye" begin
    @test bandwidths(view(Eye(∞),:,2:∞)) == (1,-1)
end

@testset "findfirst" begin
    @test findfirst(isequal(5), OneToInf()) == 5
    @test isnothing(findfirst(isequal(5.5), OneToInf()))
    @test isnothing(findfirst(isequal(-1), OneToInf()))
    @test findfirst(isequal(5), 2:∞) == 4
    @test isnothing(findfirst(isequal(5.5), 2:∞))
    @test isnothing(findfirst(isequal(-1), 2:∞))
    @test findfirst(isequal(4), 2:2:∞) == 2
    @test isnothing(findfirst(isequal(5), 2:2:∞))
    @test isnothing(findfirst(isequal(5.5), 2:2:∞))
    @test isnothing(findfirst(isequal(-1), 2:2:∞))

    @test searchsorted(Vcat(2,3:∞),10) == 9:9
    @test searchsortedfirst(Vcat(2,3:∞),10) == 9
    @test searchsortedlast(Vcat(2,3:∞),10) == 9
    @test searchsortedlast(Vcat(2,3:∞),0) == 0

    @test searchsorted(factorial.(big(1):∞), 6) == 3:3
    @test searchsortedfirst(factorial.(big(1):∞), 7) == 4
    @test searchsortedlast(factorial.(big(1):∞), 7) == 3

    @testset "Issue #178" begin
        findfirst(isone, 1:∞) == 1 
        findfirst(isone, 0:2:∞) === nothing
        findfirst(isone, -5:∞) == 7
        findfirst(isone, 2:∞) === nothing 
        findfirst(iszero, 0:∞) == 1 
        findfirst(iszero, 5:∞) === nothing 
        findfirst(iszero, 0.5:∞) === nothing 
        findfirst(iszero, -5.0:2.5:∞) == 3
    end
end

@testset "convert infrange" begin
    @test convert(AbstractArray{Float64}, 1:∞) ≡ convert(AbstractArray{Float64}, oneto(∞)) ≡
        convert(AbstractVector{Float64}, 1:∞) ≡ convert(AbstractVector{Float64}, oneto(∞)) ≡
        AbstractVector{Float64}(1:∞) ≡ AbstractVector{Float64}(oneto(∞)) ≡
        AbstractArray{Float64}(1:∞) ≡ AbstractArray{Float64}(oneto(∞)) ≡ InfUnitRange(1.0)

    @test convert(AbstractArray{Float64}, (1:∞)') ≡ convert(AbstractArray{Float64}, oneto(∞)') ≡
        convert(AbstractMatrix{Float64}, (1:∞)') ≡ convert(AbstractMatrix{Float64}, oneto(∞)') ≡
        AbstractMatrix{Float64}((1:∞)') ≡ AbstractMatrix{Float64}(oneto(∞)') ≡
        AbstractArray{Float64}((1:∞)') ≡ AbstractArray{Float64}(oneto(∞)') ≡
        InfUnitRange(1.0)'

    @test convert(AbstractArray{Float64}, transpose(1:∞)) ≡ convert(AbstractArray{Float64}, transpose(oneto(∞))) ≡
        convert(AbstractMatrix{Float64}, transpose(1:∞)) ≡ convert(AbstractMatrix{Float64}, transpose(oneto(∞))) ≡
        AbstractMatrix{Float64}(transpose(1:∞)) ≡ AbstractMatrix{Float64}(transpose(oneto(∞))) ≡
        AbstractArray{Float64}(transpose(1:∞)) ≡ AbstractArray{Float64}(transpose(oneto(∞))) ≡
        transpose(InfUnitRange(1.0))
end

@testset "cached indexing" begin
    @test cache(1:∞)[Fill(2,∞)][1:10] == fill(2,10)
    @test isassigned(cache(1:∞),RealInfinity())
    @test cache(1:∞)[RealInfinity()] == ∞
    C = cache(Fill(2.0,∞,∞))
    C[1:5,1:5] .= randn.()
    @test C[1:∞,2:∞][1:10,1:10] == C[1:10,2:11]
    @test C[1:∞,2:11][1:10,1:10] == C[1:10,2:11]
    @test C[1:10,2:∞][1:10,1:10] == C[1:10,2:11]
end

struct MyInfVector <: AbstractVector{Int} end
Base.size(::MyInfVector) = (ℵ₀,)
Base.getindex(::MyInfVector, k::Int) = k

struct MyInfMatrix <: AbstractMatrix{Int} end
Base.size(::MyInfMatrix) = (ℵ₀,ℵ₀)
Base.getindex(::MyInfMatrix, k::Int, j::Int) = k+j


@testset "MyInfArray" begin
    @test MyInfVector()[2:∞][1:10] == 2:11
    @test MyInfVector()[:][2:10] == MyInfVector()[2:10]

    @test MyInfMatrix()[2:∞,3:∞][1:10,1:10] == MyInfMatrix()[2:11,3:12]
    @test MyInfMatrix()[2:11,3:∞][1:10,1:10] == MyInfMatrix()[2:11,3:12]
    @test MyInfMatrix()[2:∞,3:12][1:10,1:10] == MyInfMatrix()[2:11,3:12]
end

struct MyReal <: Real
    x::Float64
end

Base.ArithmeticStyle(::Type{MyReal}) = Base.ArithmeticRounds()

@testset "non-float _range with ArithmeticRounds" begin
    # this missing overloaded was triggered by ForwardDiff.Dual
    @test range(MyReal(0.1); step=MyReal(0.2), length=ℵ₀) isa InfStepRange
end

@testset "UpperTriangular Inverse" begin
    A = UpperTriangular(Ones(∞,∞))
    @test ApplyArray(inv,A)[1:10,1:10] ≈ diagm(0 => ones(10), 1 => -ones(9))
end

@testset "3-mul (changed in Julia v1.9)" begin
    @test *(Eye(∞),Diagonal(1:∞), Eye(∞)) == broadcast(*, Ones(∞), Diagonal(1:∞), Ones(1,∞)) ==
        broadcast(*, Ones(∞), Diagonal(1:∞), Ones(∞,∞)) == broadcast(*, Ones(∞,1), Diagonal(1:∞), Ones(∞,∞)) ==
        broadcast(*, Ones(∞,∞), Diagonal(1:∞), Ones(∞,∞)) == Diagonal(1:∞)

    @test_throws DimensionMismatch broadcast(*, Ones(∞,2), Diagonal(1:∞))
    @test_throws DimensionMismatch broadcast(*, Diagonal(1:∞), Ones(2,∞))
end

@testset "∞-cached matrix indexing" begin
    c = zeros(5,∞)
	c[6] = 2
	@test c[6] == 2

    c = zeros(∞,3)
	c[6] = 2
	@test c[6] == 2

    c = zeros(5,3,∞);
    c[6] = 2
    @test c[1,2,1] == 2
    @test_broken c[6] == 2
end

@testset "print_matrix_row" begin
    # check that show works
    Base.sprint(show, I(ℵ₀), context=:limit=>true)
end

@testset "comprehension" begin
    @test [k * (k+1) for k = 1:∞][1:10] == [k * (k+1) for k = 1:10]
    @test [k * (k+1) for k = 1:2:∞][1:10] == [k * (k+1) for k = 1:2:20]
    @test [k * (k+1) for k = 1:2.0:∞][1:10] == [k * (k+1) for k = 1:2.0:20]
    @test [k * (k+1) for k = 1:2:∞] isa AbstractVector{Int}
    @test [k * (k+1) for k = 1:2.0:∞] isa AbstractVector{Float64}
end

@testset "diag" begin
    D = Diagonal(1:∞)
    @test @inferred(diag(D)) === 1:∞
    @test @inferred((D -> diag(D,1))(D)) === Zeros{Int}(ℵ₀)
    # test for compile-time evaluation of off-diagonals
    @inferred Val((D -> diag(D,1))(D))
    # Issue #176 
    @test inv(D)[1:100,1:100] == Diagonal(inv.(1:∞))[1:100,1:100]
end

@testset "inf padded" begin
    v = Vcat(1, Zeros(∞))
    @test LazyArrays.sub_materialize(view(v, 1:∞))[1:10] == [1; zeros(9)]
    @test LazyArrays.sub_materialize(view(v, 2:∞))[1:10] == zeros(10)
    @test v[2:∞] isa Zeros
    @test v[1:∞] == v
end

@testset "issue #180" begin
    @test isnothing(findfirst(==(21), 10:-1:-∞))
    @test isnothing(findfirst(==(11), 10:-1:-∞))
    @test findfirst(==(9), 10:-1:-∞) == 2
    r = 10:-1:-∞
    v = -20
    ind = findfirst(==(v), r)
    @test r[ind] == v
end


@testset "bounds-checking for StepRangeLen{<:CartesianIndex}" begin
    if VERSION >= v"1.11.0-rc3"
        D = Diagonal(1:∞)
        @test checkbounds(Bool, D, diagind(D, IndexCartesian()))
    end
end

@testset "Vector * ∞ matrix" begin
    a = [1+im,2+im]
    A = a * Ones{Complex{Int}}(1,∞)
    @test A[:,1:5] == a * ones(1,5)
    
    @test (a*permutedims(1:∞))[:,1:5] == a*(1:5)'
    @test (a*Hcat(Zeros(1,2), permutedims(1:∞)))[1,1:5] == (a*Vcat(Hcat(Zeros(1,2), permutedims(1:∞))))[1,1:5]
end

include("test_infconv.jl")
include("test_infblock.jl")
include("test_infbanded.jl")
include("test_infblockbanded.jl")
