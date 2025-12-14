using InfiniteArrays, BlockBandedMatrices, ArrayLayouts, FillArrays, BlockArrays, Test
using Base: oneto

@testset "BlockBanded" begin
    @testset "Triangle Recurrence" begin
        a = b = c = 0.0
        n = mortar(Fill.(oneto(∞), oneto(∞)))
        k = mortar(Base.OneTo.(oneto(∞)))
        Dy = BlockBandedMatrices._BandedBlockBandedMatrix((k .+ (b + c))', axes(k, 1), (-1, 1), (-1, 1))
        N = 100
        @test Dy[Block.(1:N), Block.(1:N)] == BlockBandedMatrices._BandedBlockBandedMatrix((k.+(b+c))[Block.(1:N)]', axes(k, 1)[Block.(1:N)], (-1, 1), (-1, 1))
        @test colsupport(Dy, axes(Dy,2)) == 1:∞
        @test rowsupport(Dy, axes(Dy,1)) == 2:∞
    end

    

    @testset "BlockTridiagonal" begin
        A = BlockTridiagonal(Vcat([fill(1.0, 2, 1), Matrix(1.0I, 2, 2), Matrix(1.0I, 2, 2), Matrix(1.0I, 2, 2)], Fill(Matrix(1.0I, 2, 2), ∞)),
            Vcat([zeros(1, 1)], Fill(zeros(2, 2), ∞)),
            Vcat([fill(1.0, 1, 2), Matrix(1.0I, 2, 2)], Fill(Matrix(1.0I, 2, 2), ∞)))

        @test isblockbanded(A)    
        @test BlockBandedMatrix(A)[1:100, 1:100] == BlockBandedMatrix(A, (2, 1))[1:100, 1:100] == BlockBandedMatrix(A, (1, 1))[1:100, 1:100] == A[1:100, 1:100]
    end
end