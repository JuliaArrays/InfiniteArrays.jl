using InfiniteArrays, DSP, FillArrays, Test
using Base: oneto

@testset "conv" begin
    @test conv(1:∞, [2]) ≡ conv([2], 1:∞) ≡ 2:2:∞
    @test conv(1:2:∞, [2]) ≡ conv([2], 1:2:∞) ≡ 2:4:∞
    @test conv(1:∞, Ones(∞))[1:5] == conv(Ones(∞),1:∞)[1:5] == [1,3,6,10,15]
    @test conv(Ones(∞), Ones(∞)) ≡ 1.0:1.0:∞
    @test conv(Ones{Int}(∞), Ones{Int}(∞)) ≡ oneto(∞)
    @test conv(Ones{Bool}(∞), Ones{Bool}(∞)) ≡ oneto(∞)
    @test conv(Fill{Int}(2,∞), Fill{Int}(1,∞)) ≡ conv(Fill{Int}(2,∞), Ones{Int}(∞)) ≡
                conv(Ones{Int}(∞), Fill{Int}(2,∞)) ≡ 2:2:∞
    @test conv(Ones{Int}(∞), [1,5,8])[1:10] == conv([1,5,8], Ones{Int}(∞))[1:10] ==
                    conv(fill(1,100),[1,5,8])[1:10] == conv([1,5,8], fill(1,100))[1:10]
    @test conv(Ones{Int}(∞), 1:4)[1:10] == conv(1:4, Ones{Int}(∞))[1:10] ==
                    conv(fill(1,100),collect(1:4))[1:10] == conv(collect(1:4), fill(1,100))[1:10]
end
