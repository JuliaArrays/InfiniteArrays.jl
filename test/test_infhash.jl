using InfiniteArrays, FillArrays, LazyArrays, ArrayLayouts, LinearAlgebra, Infinities, Test
using InfiniteArrays: InfUnitRange, InfStepRange, OneToInf, BiInfUnitRange

@testset "hash" begin
    @testset "ranges" begin
        @test hash(OneToInf()) == hash(OneToInf{Int16}()) == hash(1:∞) == hash(InfUnitRange(1)) ==
                hash(InfStepRange(1,1)) == hash(1.0:1.0:∞)
        @test hash(3:∞) == hash(InfStepRange(3,1)) == hash(3.0:∞)
        @test hash(1:2:∞) == hash(InfStepRange(1,2))
        @test hash(BiInfUnitRange()) == hash(BiInfUnitRange{Float64}())
        @test hash(view(1:∞, 2:∞)) == hash(2:∞)

        # `Base` axis wrappers compare equal to the range they wrap, so they have to hash alike
        @test hash(Base.Slice(OneToInf())) == hash(Base.IdentityUnitRange(1:∞)) == hash(1:∞)

        rs = (1:∞, 2:∞, 1:2:∞, 2:2:∞, 1.5:1.0:∞, BiInfUnitRange())
        @test length(unique(hash.(rs))) == length(rs)
        for r in rs, h in (zero(UInt), UInt(5))
            @test hash(r, h) isa UInt
        end
    end

    @testset "wrappers" begin
        # `Diagonal(1:∞) == Diagonal(1.0:1.0:∞)`, so the two have to hash alike
        @test hash(Diagonal(1:∞)) == hash(Diagonal(1.0:1.0:∞))
        @test hash(Diagonal(1:∞)) ≠ hash(Diagonal(2:∞))
        @test hash(Diagonal(1:∞)') == hash(Diagonal(1:∞))
        @test hash((1:∞)') == hash(transpose(1:∞))
        @test hash((1:∞)') ≠ hash(1:∞)
        @test hash(Tridiagonal(1:∞, 1:∞, 1:∞)) isa UInt
        @test hash(SymTridiagonal(1:∞, 1:∞)) isa UInt
        @test hash(Bidiagonal(1:∞, 1:∞, :U)) isa UInt
        @test hash(Symmetric(Diagonal(1:∞))) isa UInt
        @test hash(UpperTriangular(Diagonal(1:∞))) isa UInt
        for h in (zero(UInt), UInt(5))
            @test hash(Diagonal(1:∞), h) == hash(Diagonal(1.0:1.0:∞), h)
        end

        # a view of an infinite range need not be infinite
        @test hash(view(1:∞, 1:5)) == hash(1:5)
        @test hash(view(1:∞, 1:5), UInt(7)) == hash(1:5, UInt(7))
    end

    @testset "isequal" begin
        @test isequal(1:∞, OneToInf())
        @test isequal(1:∞, InfStepRange(1,1))
        @test isequal(1:∞, 1.0:1.0:∞)
        @test isequal(Base.Slice(OneToInf()), 1:∞)
        @test isequal(BiInfUnitRange(), BiInfUnitRange{Float64}())
        @test !isequal(1:∞, 2:∞)
        @test !isequal(1:∞, 1:2:∞)
        @test !isequal(1:∞, BiInfUnitRange())
        # offset wrappers hold the same entries but are not `isequal`, as in `Base`
        @test !isequal(Base.IdentityUnitRange(3:∞), 3:∞)
        @test hash(Base.IdentityUnitRange(3:∞)) ≠ hash(3:∞)
        @test isequal(Base.IdentityUnitRange(3:∞), Base.IdentityUnitRange(3:∞))

        # `isequal` and `hash` together are what `Dict` and `Set` need
        d = Dict(1:∞ => "a", 1:2:∞ => "b", 3:∞ => "c")
        @test d[OneToInf()] == "a"
        @test d[InfStepRange(1,2)] == "b"
        @test d[InfUnitRange(3)] == "c"
        @test length(Set([1:∞, OneToInf(), InfUnitRange(1), InfStepRange(1,1)])) == 1
    end

    @testset "finite ranges hash as in Base" begin
        # nothing here should divert a finite range from the `Base` implementation
        for (a, b) in ((Base.Slice(Base.OneTo(3)), 1:3),
                       (Base.IdentityUnitRange(Base.OneTo(4)), 1:4))
            @test hash(a) == hash(b)
            @test hash(a, UInt(7)) == hash(b, UInt(7))
        end
    end
end
