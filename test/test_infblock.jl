using InfiniteArrays, BlockArrays, LazyArrays, Test
using InfiniteArrays: OneToInf, RealInfinity
using Base: oneto
using LazyArrays: LazyArrayStyle
using BlockArrays: BlockTridiagonal

const InfiniteArraysBlockArraysExt = Base.get_extension(InfiniteArrays, :InfiniteArraysBlockArraysExt)




@testset "∞-block arrays" begin
    @testset "blockedonetoinf" begin
        b = blockedrange(OneToInf())
        b2 = b .+ b;
        for i in 1:10
            @test b2[Block(i)] == b[Block(i)] + b[Block(i)]
        end
    end

    @testset "fixed block size" begin
        k = Base.OneTo.(oneto(∞))
        n = Fill.(oneto(∞), oneto(∞))
        @test broadcast(length, k) ≡ map(length, k) ≡ OneToInf()
        @test broadcast(length, n) ≡ map(length, n) ≡ OneToInf()

        b = mortar(Fill([1, 2], ∞))
        @test blockaxes(b, 1) ≡ Block.(OneToInf())
        @test b[Block(5)] == [1, 2]
        @test b[Block.(2:∞)][Block.(2:10)] == b[Block.(3:11)]
        @test exp.(b)[Block.(2:∞)][Block.(2:10)] == exp.(b[Block.(3:11)])

        @test blockedrange(Vcat(2, Fill(3, ∞))) isa BlockedOneTo{<:Any,<:InfiniteArrays.InfStepRange}

        c = BlockedArray(1:∞, Vcat(2, Fill(3, ∞)))
        @test c[Block.(2:∞)][Block.(2:10)] == c[Block.(3:11)]

        @test length(axes(b, 1)) ≡ ℵ₀
        @test last(axes(b, 1)) ≡ ℵ₀
        @test Base.BroadcastStyle(typeof(b)) isa LazyArrayStyle{1}
    end

    @testset "1:∞ blocks" begin
        a = blockedrange(oneto(∞))
        @test axes(a, 1) == a
        o = Ones((a,))
        @test Base.BroadcastStyle(typeof(a)) isa LazyArrayStyle{1}
        b = exp.(a)
        @test axes(b, 1) == a
        @test o .* b isa typeof(b)
        @test b .* o isa typeof(b)
    end

    @testset "padded" begin
        c = BlockedArray([1; zeros(∞)], Vcat(2, Fill(3, ∞)))
        @test c + c isa BlockedVector
    end


    @testset "triangle recurrences" begin
        @testset "n and k" begin
            n = mortar(Fill.(oneto(∞), oneto(∞)))
            k = mortar(Base.OneTo.(oneto(∞)))

            @test n[Block(5)] ≡ layout_getindex(n, Block(5)) ≡ view(n, Block(5)) ≡ Fill(5, 5)
            @test k[Block(5)] ≡ layout_getindex(k, Block(5)) ≡ view(k, Block(5)) ≡ Base.OneTo(5)
            @test Base.BroadcastStyle(typeof(n)) isa LazyArrays.LazyArrayStyle{1}
            @test Base.BroadcastStyle(typeof(k)) isa LazyArrays.LazyArrayStyle{1}

            N = 1000
            v = view(n, Block.(Base.OneTo(N)))
            @test view(v, Block(2)) ≡ Fill(2, 2)
            @test axes(v) isa Tuple{BlockedOneTo{Int,ArrayLayouts.RangeCumsum{Int64,Base.OneTo{Int64}}}}
            @test @allocated(axes(v)) ≤ 40

            dest = BlockedArray{Float64}(undef, axes(v))
            @test copyto!(dest, v) == v
            @test @allocated(copyto!(dest, v)) ≤ 40

            v = view(k, Block.(Base.OneTo(N)))
            @test view(v, Block(2)) ≡ Base.OneTo(2)
            @test axes(v) isa Tuple{BlockedOneTo{Int,ArrayLayouts.RangeCumsum{Int64,Base.OneTo{Int64}}}}
            @test @allocated(axes(v)) ≤ 40
            @test copyto!(dest, v) == v

            @testset "stack overflow" begin
                i = Base.to_indices(k, (Block.(2:∞),))[1].indices
                @test last(i) == ℵ₀
            end

            v = view(k, Block.(2:∞))
            @test Base.BroadcastStyle(typeof(v)) isa LazyArrayStyle{1}
            @test v[Block(1)] == 1:2
            @test v[Block(1)] ≡ k[Block(2)] ≡ Base.OneTo(2)

            @test axes(n, 1) isa BlockedOneTo{Int,ArrayLayouts.RangeCumsum{Int64,OneToInf{Int64}}}
        end
    end

    @testset "blockdiag" begin
        D = Diagonal(mortar(Fill.((-(0:∞) - (0:∞) .^ 2), 1:2:∞)))
        x = [randn(5); zeros(∞)]
        x̃ = BlockedArray(x, (axes(D, 1),))
        @test (D*x)[1:10] == (D*x̃)[1:10]
    end

    @testset "sortedunion" begin
        a = cumsum(1:2:∞)
        @test BlockArrays.sortedunion(a, a) ≡ a
        @test BlockArrays.sortedunion([∞], a) ≡ BlockArrays.sortedunion(a, [∞]) ≡ a
        @test BlockArrays.sortedunion([∞], [∞]) == [∞]

        b = Vcat([1, 2], 3:∞)
        c = Vcat(1, 3:∞)
        @test BlockArrays.sortedunion(b, b) ≡ b
        @test BlockArrays.sortedunion(c, c) ≡ c
    end

    @testset "Algebra" begin
        @testset "Triangle OP recurrences" begin
            k = mortar(Base.OneTo.(1:∞))
            n = mortar(Fill.(1:∞, 1:∞))
            @test k[Block.(2:3)] isa BlockArray
            @test n[Block.(2:3)] isa BlockArray
            @test k[Block.(2:3)] == [1, 2, 1, 2, 3]
            @test n[Block.(2:3)] == [2, 2, 3, 3, 3]
            @test blocksize(BroadcastVector(exp, k)) == (ℵ₀,)
            @test BroadcastVector(exp, k)[Block.(2:3)] == exp.([1, 2, 1, 2, 3])
            # BroadcastVector(+,k,n)
        end
        # Multivariate OPs Corollary (3)
        # n = 5
        # BlockTridiagonal(Zeros.(1:∞,2:∞),
        #         (n -> Diagonal(((n+2).+(0:n)))/ (2n + 2)).(0:∞),
        #         Zeros.(2:∞,1:∞))
    end
    
    @testset "findblock at +∞, HarmonicOrthogonalPolynomials#88" begin
        @test findblock(blockedrange(1:2:∞), RealInfinity()) == Block(ℵ₀)
    end    
end