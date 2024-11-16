using InfiniteArrays, BlockArrays
using InfiniteArrays: OneToInf

@testset "blockedonetoinf" begin
    b = blockedrange(OneToInf())
    b2 = b .+ b;
    for i in 1:10
        @test b2[Block(i)] == b[Block(i)] + b[Block(i)]
    end
end