using ArrayLayouts, InfiniteArrays, BandedMatrices, FillArrays, LazyArrays, Test
import BandedMatrices: _BandedMatrix, bandeddata

const InfiniteArraysBandedMatricesExt = Base.get_extension(InfiniteArrays, :InfiniteArraysBandedMatricesExt)
const InfToeplitz = InfiniteArraysBandedMatricesExt.InfToeplitz
const TriToeplitz = InfiniteArraysBandedMatricesExt.TriToeplitz
const SymTriPertToeplitz = InfiniteArraysBandedMatricesExt.SymTriPertToeplitz
const TriPertToeplitz = InfiniteArraysBandedMatricesExt.TriPertToeplitz
const AdjTriPertToeplitz = InfiniteArraysBandedMatricesExt.AdjTriPertToeplitz
const ConstRows = InfiniteArraysBandedMatricesExt.ConstRows
const BandedToeplitzLayout = InfiniteArraysBandedMatricesExt.BandedToeplitzLayout
const TridiagonalToeplitzLayout = InfiniteArraysBandedMatricesExt.TridiagonalToeplitzLayout
const BidiagonalToeplitzLayout = InfiniteArraysBandedMatricesExt.BidiagonalToeplitzLayout
const PertToeplitz = InfiniteArraysBandedMatricesExt.PertToeplitz
const PertToeplitzLayout = InfiniteArraysBandedMatricesExt.PertToeplitzLayout
const InfBandCartesianIndices = InfiniteArraysBandedMatricesExt.InfBandCartesianIndices

using Base: oneto
using LazyArrays: simplifiable, ApplyLayout, BroadcastBandedLayout

@testset "∞-banded" begin
    @testset "Diagonal and BandedMatrix" begin
        D = Diagonal(Fill(2,∞))

        B = D[1:∞,2:∞]
        @test B isa BandedMatrix
        @test B[1:10,1:10] == diagm(-1 => Fill(2,9))
        @test B[1:∞,2:∞] isa BandedMatrix

        A = BandedMatrix(0 => 1:∞, 1=> Fill(2.0,∞), -1 => Fill(3.0,∞))
        x = [1; 2; zeros(∞)]
        @test A*x isa Vcat
        @test (A*x)[1:10] == A[1:10,1:10]*x[1:10]

        @test InfBandCartesianIndices(0)[1:5] == CartesianIndex.(1:5,1:5)
        @test InfBandCartesianIndices(1)[1:5] == CartesianIndex.(1:5,2:6)
        @test InfBandCartesianIndices(-1)[1:5] == CartesianIndex.(2:6,1:5)

        @test A[band(0)][2:10] == 2:10
        @test D[band(0)] ≡ Fill(2,∞)
        @test D[band(1)] ≡ Fill(0,∞)

        @test B*A*x isa Vcat
        @test (B*A*x)[1:10] == [0; 10; 14; 12; zeros(6)]

        @test _BandedMatrix((1:∞)', ∞, -1, 1) isa BandedMatrix
    end

    @testset "∞-Toeplitz" begin
        A = BandedMatrix(1 => Fill(2im,∞), 2 => Fill(-1,∞), 3 => Fill(2,∞), -2 => Fill(-4,∞), -3 => Fill(-2im,∞))
        @test A isa InfToeplitz
        @test MemoryLayout(A.data) == ConstRows()
        @test MemoryLayout(A) == BandedToeplitzLayout()
        @test LazyArrays.islazy(A) == Val(true)

        V = view(A,:,3:∞)
        @test MemoryLayout(typeof(bandeddata(V))) == ConstRows()
        @test MemoryLayout(typeof(V)) == BandedToeplitzLayout()
        @test BandedMatrix(V) isa InfToeplitz
        @test A[:,3:end] isa InfToeplitz

        @test (A + 2I)[1:10,1:10] == (2I + A)[1:10,1:10] == A[1:10,1:10] + 2I
        @test (A*A)[1:10,1:10] ≈ A[1:10,1:16]A[1:16,1:10]
        @test (A * Fill(2,∞))[1:10] ≈ 2A[1:10,1:16]*ones(16)
        @test (Fill(2,∞,∞)*A)[1:10,1:10] ≈ fill(2,10,13)A[1:13,1:10]

        @test Eye(∞) * A isa BandedMatrix
        @test A * Eye(∞) isa BandedMatrix

        @test A * [1; 2; Zeros(∞)] isa Vcat
        @test A * [1; 2; Zeros(∞)] == [A[1:5,1:2] * [1,2]; Zeros(∞)]

        @test A * Vcat([1 2; 3 4], Zeros(∞,2)) isa ApplyMatrix{ComplexF64,typeof(Base.setindex)}
        @test simplifiable(*, A, Vcat([1 2; 3 4], Zeros(∞,2))) == Val(true)

        @test MemoryLayout(Tridiagonal(Fill(1,∞), Fill(2,∞), Fill(3,∞))) isa TridiagonalToeplitzLayout
        @test MemoryLayout(Bidiagonal(Fill(1,∞), Fill(2,∞), :U)) isa BidiagonalToeplitzLayout
        @test MemoryLayout(SymTridiagonal(Fill(1,∞), Fill(2,∞))) isa TridiagonalToeplitzLayout


        @testset "algebra" begin
            T = Tridiagonal(Fill(1,∞), Fill(2,∞), Fill(3,∞))
            @test T isa TriToeplitz
            @test (T + 2I)[1:10,1:10] == (2I + T)[1:10,1:10] == T[1:10,1:10] + 2I
        end

        @testset "constant data" begin
            A = BandedMatrix(1 => Fill(2im,∞), 2 => Fill(-1,∞), 3 => Fill(2,∞), -2 => Fill(-4,∞), -3 => Fill(-2im,∞))
            B = _BandedMatrix(Fill(2,4,∞), ∞, 1,2)
            @test (B*B)[1:10,1:10] ≈ B[1:10,1:14]B[1:14,1:10]
            @test (A*B)[1:10,1:10] ≈ A[1:10,1:14]B[1:14,1:10]
            @test (B*A)[1:10,1:10] ≈ B[1:10,1:14]A[1:14,1:10]
            @test simplifiable(*, B, B) == Val(true)
            @test length((A*B*B).args) == 2
            @test length((B*B*A).args) == 2
        end

        @testset "Toep * Diag" begin
            A = BandedMatrix(1 => Fill(2im,∞), 2 => Fill(-1,∞), 3 => Fill(2,∞), -2 => Fill(-4,∞), -3 => Fill(-2im,∞))
            D = Diagonal(1:∞)
            @test D*A isa BroadcastMatrix
            @test A*D isa BroadcastMatrix
            @test simplifiable(*, D, A) == Val(true)
            @test simplifiable(*, A, D) == Val(true)
        end
    end

    @testset "Pert-Toeplitz" begin
        @testset "Inf Pert" begin
            A = BandedMatrix(-2 => Vcat(Float64[], Fill(1/4,∞)), 0 => Vcat([1.0+im,2,3],Fill(0,∞)), 1 => Vcat(Float64[], Fill(1,∞)))
            @test A isa PertToeplitz
            @test MemoryLayout(A) isa PertToeplitzLayout
            V = view(A,2:∞,2:∞)
            @test MemoryLayout(V) isa PertToeplitzLayout
            @test BandedMatrix(V) isa PertToeplitz
            @test A[2:∞,2:∞] isa PertToeplitz

            @test (A + 2I)[1:10,1:10] == (2I + A)[1:10,1:10] == A[1:10,1:10] + 2I

            @test Eye(∞) * A isa BandedMatrix
            @test A * Eye(∞) isa BandedMatrix
        end

        @testset "TriPert" begin
            A = SymTridiagonal(Vcat([1,2.], Fill(2.,∞)), Vcat([3.,4.], Fill.(0.5,∞)))
            @test A isa SymTriPertToeplitz
            @test (A + 2I)[1:10,1:10] == (2I + A)[1:10,1:10] == A[1:10,1:10] + 2I

            A = Tridiagonal(Vcat([3.,4.], Fill.(0.5,∞)), Vcat([1,2.], Fill(2.,∞)), Vcat([3.,4.], Fill.(0.5,∞)))
            @test A isa TriPertToeplitz
            @test Adjoint(A) isa AdjTriPertToeplitz
            @test (A + 2I)[1:10,1:10] == (2I + A)[1:10,1:10] == A[1:10,1:10] + 2I
            @test (Adjoint(A) + 2I)[1:10,1:10] == (2I + Adjoint(A))[1:10,1:10] == Adjoint(A)[1:10,1:10] + 2I
        end


        @testset "InfBanded" begin
            A = _BandedMatrix(Fill(2,4,∞),ℵ₀,2,1)
            B = _BandedMatrix(Fill(3,2,∞),ℵ₀,-1,2)
            @test mul(A,A) isa PertToeplitz
            @test A*A isa PertToeplitz
            @test (A*A)[1:20,1:20] == A[1:20,1:23]*A[1:23,1:20]
            @test (A*B)[1:20,1:20] == A[1:20,1:23]*B[1:23,1:20]
        end
    end

    @testset "adjortrans" begin
        A = BandedMatrix(0 => 1:∞, 1=> Fill(2.0+im,∞), -1 => Fill(3.0,∞))
        @test copy(A')[1:10,1:10] == (A')[1:10,1:10]
        @test copy(transpose(A))[1:10,1:10] == transpose(A)[1:10,1:10]
    end

    @testset "Eye subindex" begin
        @test Eye(∞)[:,1:3][1:5,:] == Eye(∞)[Base.Slice(oneto(∞)),1:3][1:5,:] == Eye(5,3)
        @test Eye(∞)[1:3,:][:,1:5] == Eye(∞)[1:3,Base.Slice(oneto(∞))][:,1:5] == Eye(3,5)
        @test Eye(∞)[:,:][1:5,1:3] == Eye(∞)[Base.Slice(oneto(∞)),Base.Slice(oneto(∞))][1:5,1:3] == Eye(5,3)
    end

    @testset "band(0) indexing" begin
        D = ApplyArray(*, Diagonal(1:∞), Diagonal(1:∞))
        @test D[band(0)][1:10] == (1:10).^2
    end

    @testset "Fill * Banded" begin
        A = _BandedMatrix(Ones(1,∞), ∞, 1,-1)
        B = _BandedMatrix(Fill(1.0π,1,∞), ∞, 0,0)
        @test (A*B)[1:10,1:10] ≈ BandedMatrix(-1 => Fill(1.0π,9))
    end

    @testset "concat" begin
        H = ApplyArray(hvcat, 2, 1, [1 Zeros(1,∞)], [1; Zeros(∞)], Diagonal(1:∞))
        @test bandwidths(H) == (1,1)
        H = ApplyArray(hvcat, 2, 1, [0 Zeros(1,∞)], [0; Zeros(∞)], Diagonal(1:∞))
        @test bandwidths(H) == (0,0)
        H = ApplyArray(hvcat, (2,2), 1, [1 Zeros(1,∞)], [1; Zeros(∞)], Diagonal(1:∞))
        @test_broken bandwidths(H) == (1,1)
    end

    @testset "Banded * PaddedMatrix" begin
        A = Eye(∞)[2:∞,:]
        @test LazyArrays.islazy(A) == Val(true)
        B = PaddedArray(randn(3,3),ℵ₀,ℵ₀)
        @test (A*B)[1:10,1:10] ≈ A[1:10,1:10] * B[1:10,1:10]
    end

    @testset "SubArray broadcasting" begin
        A = BandedMatrix(2 => 1:∞)
        @test exp.(A[1:2:∞,1:2:∞])[1:10,1:10] ≈ exp.(A[1:2:20,1:2:20])
        @test A[band(2)][1:5] == 1:5
        @test _BandedMatrix((1:∞)', ∞, -1,1)[band(1)][1:5] == 2:6
        @test exp.(view(A,band(2)))[1:10] ≈ exp.(1:10)

        @test BandedMatrices.banded_similar(Int, (∞,5), (1,1)) isa BandedMatrix
        @test BandedMatrices.banded_similar(Int, (5,∞), (1,1)) isa Adjoint{<:Any,<:BandedMatrix}

        A = BandedMatrix{Int}((2 => 1:∞,), (∞,∞), (0,2))
        @test eltype(A) == Int
        @test bandwidths(A) == (0,2)

        A = BandedMatrix{Int}((2 => Vcat([1,2], Fill(2,∞)),), (∞,∞), (0,2))
        @test A[band(2)][1:5] == [1; fill(2,4)]
    end

    @testset "Algebra" begin
        A = BandedMatrix(-3 => Fill(7 / 10, ∞), -2 => 1:∞, 1 => Fill(2im, ∞))
        @test A isa BandedMatrix{ComplexF64}
        @test A[1:10, 1:10] == diagm(-3 => Fill(7 / 10, 7), -2 => 1:8, 1 => Fill(2im, 9))

        A = BandedMatrix(0 => Vcat([1, 2, 3], Zeros(∞)), 1 => Vcat(1, Zeros(∞)))
        @test A[1, 2] == 1

        A = BandedMatrix(-3 => Fill(7 / 10, ∞), -2 => Fill(1, ∞), 1 => Fill(2im, ∞))
        Ac = BandedMatrix(A')
        At = BandedMatrix(transpose(A))
        @test Ac[1:10, 1:10] ≈ (A')[1:10, 1:10] ≈ A[1:10, 1:10]'
        @test At[1:10, 1:10] ≈ transpose(A)[1:10, 1:10] ≈ transpose(A[1:10, 1:10])

        A = BandedMatrix(-1 => Vcat(Float64[], Fill(1 / 4, ∞)), 0 => Vcat([1.0 + im], Fill(0, ∞)), 1 => Vcat(Float64[], Fill(1, ∞)))
        @test MemoryLayout(typeof(view(A.data, :, 1:10))) == ApplyLayout{typeof(hcat)}()
        Ac = BandedMatrix(A')
        At = BandedMatrix(transpose(A))
        @test Ac[1:10, 1:10] ≈ (A')[1:10, 1:10] ≈ A[1:10, 1:10]'
        @test At[1:10, 1:10] ≈ transpose(A)[1:10, 1:10] ≈ transpose(A[1:10, 1:10])

        A = BandedMatrix(-2 => Vcat(Float64[], Fill(1 / 4, ∞)), 0 => Vcat([1.0 + im, 2, 3], Fill(0, ∞)), 1 => Vcat(Float64[], Fill(1, ∞)))
        Ac = BandedMatrix(A')
        At = BandedMatrix(transpose(A))
        @test Ac[1:10, 1:10] ≈ (A')[1:10, 1:10] ≈ A[1:10, 1:10]'
        @test At[1:10, 1:10] ≈ transpose(A)[1:10, 1:10] ≈ transpose(A[1:10, 1:10])

        A = _BandedMatrix(Fill(1, 4, ∞), ℵ₀, 1, 2)
        @test A^2 isa BandedMatrix
        @test (A^2)[1:10, 1:10] == (A*A)[1:10, 1:10] == (A[1:100, 1:100]^2)[1:10, 1:10]
        @test A^3 isa ApplyMatrix{<:Any,typeof(*)}
        @test (A^3)[1:10, 1:10] == (A*A*A)[1:10, 1:10] == ((A*A)*A)[1:10, 1:10] == (A*(A*A))[1:10, 1:10] == (A[1:100, 1:100]^3)[1:10, 1:10]

        @testset "∞ x finite" begin
            A = BandedMatrix(1 => 1:∞) + BandedMatrix(-1 => Fill(2, ∞))
            B = _BandedMatrix(randn(3, 5), ℵ₀, 1, 1)

            @test lmul!(2.0, copy(B)')[:, 1:10] == (2B')[:, 1:10]

            @test_throws ArgumentError BandedMatrix(A)
            @test A * B isa MulMatrix
            @test B'A isa MulMatrix

            @test all(diag(A[1:6, 1:6]) .=== zeros(Int,6))

            @test (A*B)[1:7, 1:5] ≈ A[1:7, 1:6] * B[1:6, 1:5]
            @test (B'A)[1:5, 1:7] ≈ (B')[1:5, 1:6] * A[1:6, 1:7]
        end
    end

    @testset "Fill" begin
        A = _BandedMatrix(Ones(1, ∞), ℵ₀, -1, 1)
        @test 1.0 .* A isa BandedMatrix{Float64,<:Fill}
        @test Zeros(∞) .* A ≡ Zeros(∞, ∞) .* A ≡ A .* Zeros(1, ∞) ≡ A .* Zeros(∞, ∞) ≡ Zeros(∞, ∞)
        @test Ones(∞) .* A isa BandedMatrix{Float64,<:Ones}
        @test A .* Ones(1, ∞) isa BandedMatrix{Float64,<:Ones}
        @test 2.0 .* A isa BandedMatrix{Float64,<:Fill}
        @test A .* 2.0 isa BandedMatrix{Float64,<:Fill}
        @test Eye(∞) * A isa BandedMatrix{Float64,<:Ones}
        @test A * Eye(∞) isa BandedMatrix{Float64,<:Ones}

        @test A * A isa BandedMatrix
        @test (A*A)[1:10, 1:10] == BandedMatrix(2 => Ones(8))

        Ã = _BandedMatrix(Fill(1, 1, ∞), ℵ₀, -1, 1)
        @test A * Ã isa BandedMatrix
        @test Ã * A isa BandedMatrix
        @test Ã * Ã isa BandedMatrix

        B = _BandedMatrix(Ones(1, 10), ℵ₀, -1, 1)
        C = _BandedMatrix(Ones(1, 10), 10, -1, 1)
        D = _BandedMatrix(Ones(1, ∞), 10, -1, 1)

        @test (A*B)[1:10, 1:10] == (B*C)[1:10, 1:10] == (D*A)[1:10, 1:10] == D * B == (C*D)[1:10, 1:10] == BandedMatrix(2 => Ones(8))
    end

    @testset "Banded Broadcast" begin
        A = _BandedMatrix((1:∞)', ℵ₀, -1, 1)
        @test 2.0 .* A isa BandedMatrix{Float64,<:Adjoint}
        @test A .* 2.0 isa BandedMatrix{Float64,<:Adjoint}
        @test Eye(∞) * A isa BandedMatrix{Float64,<:Adjoint}
        @test A * Eye(∞) isa BandedMatrix{Float64,<:Adjoint}
        A = _BandedMatrix(Vcat((1:∞)', Ones(1, ∞)), ℵ₀, 0, 1)
        @test 2.0 .* A isa BandedMatrix
        @test A .* 2.0 isa BandedMatrix
        @test Eye(∞) * A isa BandedMatrix
        @test A * Eye(∞) isa BandedMatrix
        b = 1:∞
        @test bandwidths(b .* A) == (0, 1)

        @test colsupport(b .* A, 1) == 1:1
        @test Base.replace_in_print_matrix(b .* A, 2, 1, "0.0") == " ⋅ "
        @test bandwidths(A .* b) == (0, 1)
        @test A .* b' isa BroadcastArray
        @test bandwidths(A .* b') == bandwidths(A .* b')
        @test colsupport(A .* b', 3) == 2:3

        A = _BandedMatrix(Ones{Int}(1, ∞), ℵ₀, 0, 0)'
        B = _BandedMatrix((-2:-2:-∞)', ℵ₀, -1, 1)
        C = Diagonal(2 ./ (1:2:∞))
        @test bandwidths(A * (B * C)) == (-1, 1)
        @test bandwidths((A * B) * C) == (-1, 1)

        A = _BandedMatrix(Ones{Int}(1, ∞), ℵ₀, 0, 0)'
        B = _BandedMatrix((-2:-2:-∞)', ℵ₀, -1, 1)
        @test MemoryLayout(A + B) isa BroadcastBandedLayout{typeof(+)}
        @test MemoryLayout(2 * (A + B)) isa BroadcastBandedLayout{typeof(*)}
        @test bandwidths(A + B) == (0, 1)
        @test bandwidths(2 * (A + B)) == (0, 1)
    end

end
