using InfiniteArrays, BandedMatrices, Base64, Test
using InfiniteArrays: BiInfUnitRange

@testset "-∞:∞" begin
    r = BiInfUnitRange()
    @test stringmime("text/plain", r) == "BiInfUnitRange()"
    @test Fill(1, (r,))[-5] == 1
    @test exp.(r)[-3] == exp(-3)
    @test stringmime("text/plain", r.^2; context= IOContext(IOBuffer(), :limit=>true)) == "(ℵ₀-element BiInfUnitRange{$Int} with indices BiInfUnitRange()) .^ 2 with indices BiInfUnitRange():\n ⋮\n 25\n 16\n  9\n  4\n  1\n \e[1m 0\e[0m\n  1\n  4\n  9\n 16\n 25\n  ⋮"

    @test isassigned(r', 1,-12)
    @test stringmime("text/plain", r'; context= IOContext(IOBuffer(), :limit=>true)) == "1×ℵ₀ adjoint(::BiInfUnitRange{$Int}) with eltype Int64 with indices Base.OneTo(1)×BiInfUnitRange():\n  …   -6  -5  -4  -3  -2  -1  0  1  2  3  …  "

    # _BandedMatrix(
end