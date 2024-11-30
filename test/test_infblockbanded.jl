const BlockTriPertToeplitz = InfiniteArraysBlockArraysExt.BlockTriPertToeplitz

@testset "BlockBanded" begin
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
    
            @test A isa BlockTriPertToeplitz
            @test isblockbanded(A)
    
            @test A[Block.(1:2), Block(1)] == A[1:3, 1:1] == reshape([0.0, 1.0, 1.0], 3, 1)
    
            @test BlockBandedMatrix(A)[1:100, 1:100] == BlockBandedMatrix(A, (2, 1))[1:100, 1:100] == BlockBandedMatrix(A, (1, 1))[1:100, 1:100] == A[1:100, 1:100]
    
            @test (A-I)[1:100, 1:100] == A[1:100, 1:100] - I
            @test (A+I)[1:100, 1:100] == A[1:100, 1:100] + I
            @test (I+A)[1:100, 1:100] == I + A[1:100, 1:100]
            @test (I-A)[1:100, 1:100] == I - A[1:100, 1:100]
    
            @test (A-im*I)[1:100, 1:100] == A[1:100, 1:100] - im * I
            @test (A+im*I)[1:100, 1:100] == A[1:100, 1:100] + im * I
            @test (im*I+A)[1:100, 1:100] == im * I + A[1:100, 1:100]
            @test (im*I-A)[1:100, 1:100] == im * I - A[1:100, 1:100]
    
            T = mortar(LazyBandedMatrices.Tridiagonal(Fill([1 2; 3 4], ∞), Fill([1 2; 3 4], ∞), Fill([1 2; 3 4], ∞)));
            #TODO: copy BlockBidiagonal code from BlockBandedMatrices to LazyBandedMatrices
            @test T[Block(2, 2)] == [1 2; 3 4]
            @test_broken T[Block(1, 3)] == Zeros(2, 2)
        end