using InfiniteArrays, Base64, Test
using InfiniteArrays: BiInfUnitRange

@testset "-∞:∞" begin
    r = BiInfUnitRange()
    @test stringmime("text/plain", r) == "BiInfUnitRange()"
    @test Fill(1, (r,))[-5] == 1
    @test exp.(r)[-3] == exp(-3)
end