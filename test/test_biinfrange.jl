using InfiniteArrays, Base64, Test
using InfiniteArrays: BiInfUnitRange

@testset "-∞:∞" begin
    r = BiInfUnitRange()
    
    @test AbstractArray{Int}(r) ≡ AbstractVector{Int}(r) ≡ r
    @test Base.has_offset_axes(r)

    @test r[-∞] ≡ -∞
    @test r[∞] ≡ r[ℵ₀] ≡ ℵ₀
    @test r[-3:3] ≡ -3:3
    
    
    
    @test stringmime("text/plain", r) == "BiInfUnitRange()"
    
    @test Fill(1, (r,))[-5] == 1
    @test exp.(r)[-3] == exp(-3)
    
    
    @test stringmime("text/plain", r * (1:10)'; context= IOContext(IOBuffer(), :limit=>true)) == "(ℵ₀-element BiInfUnitRange{$Int} with indices BiInfUnitRange()) .* (1×10 adjoint(::UnitRange{$Int}) with eltype $Int) with indices BiInfUnitRange()×Base.OneTo(10):\n ⋮                        ⋮                 \n -5  -10  -15  -20  -25  -30  -35  -40  -45  -50\n -4   -8  -12  -16  -20  -24  -28  -32  -36  -40\n -3   -6   -9  -12  -15  -18  -21  -24  -27  -30\n -2   -4   -6   -8  -10  -12  -14  -16  -18  -20\n -1   -2   -3   -4   -5   -6   -7   -8   -9  -10\n \e[1m 0    0    0    0    0    0    0    0    0    0\e[0m\n  1    2    3    4    5    6    7    8    9   10\n  2    4    6    8   10   12   14   16   18   20\n  3    6    9   12   15   18   21   24   27   30\n  4    8   12   16   20   24   28   32   36   40\n  5   10   15   20   25   30   35   40   45   50\n  ⋮                        ⋮                 "
    @test stringmime("text/plain", r * (1:1000)'; context= IOContext(IOBuffer(), :limit=>true)) == "(ℵ₀-element BiInfUnitRange{$Int} with indices BiInfUnitRange()) .* (1×1000 adjoint(::UnitRange{$Int}) with eltype $Int) with indices BiInfUnitRange()×Base.OneTo(1000):
 -5  -10  -15  -20  -25  -30  -35  -40  …  -4980  -4985  -4990  -4995  -5000
 -4   -8  -12  -16  -20  -24  -28  -32     -3984  -3988  -3992  -3996  -4000
 -3   -6   -9  -12  -15  -18  -21  -24     -2988  -2991  -2994  -2997  -3000
 -2   -4   -6   -8  -10  -12  -14  -16     -1992  -1994  -1996  -1998  -2000
 -1   -2   -3   -4   -5   -6   -7   -8      -996   -997   -998   -999  -1000
  0    0    0    0    0    0    0    0  …      0      0      0      0      0
  1    2    3    4    5    6    7    8       996    997    998    999   1000
  2    4    6    8   10   12   14   16      1992   1994   1996   1998   2000
  3    6    9   12   15   18   21   24      2988   2991   2994   2997   3000
  4    8   12   16   20   24   28   32      3984   3988   3992   3996   4000
  5   10   15   20   25   30   35   40  …   4980   4985   4990   4995   5000
  ⋮                        ⋮            ⋱      ⋮                       "
    
    @test stringmime("text/plain", r.^2; context= IOContext(IOBuffer(), :limit=>true)) == "(ℵ₀-element BiInfUnitRange{$Int} with indices BiInfUnitRange()) .^ 2 with indices BiInfUnitRange():\n ⋮\n 25\n 16\n  9\n  4\n  1\n \e[1m 0\e[0m\n  1\n  4\n  9\n 16\n 25\n  ⋮"
    
    @test isassigned(r', 1,-12)
    @test stringmime("text/plain", r'; context= IOContext(IOBuffer(), :limit=>true)) == "1×ℵ₀ adjoint(::BiInfUnitRange{$Int}) with eltype $Int with indices Base.OneTo(1)×BiInfUnitRange():\n  …   -6  -5  -4  -3  -2  -1  0  1  2  3  …  "
    
    @test stringmime("text/plain", r * r'; context= IOContext(IOBuffer(), :limit=>true)) == "(ℵ₀-element BiInfUnitRange{$Int} with indices BiInfUnitRange()) .* (1×ℵ₀ adjoint(::BiInfUnitRange{$Int}) with eltype $Int with indices Base.OneTo(1)×BiInfUnitRange()) with indices BiInfUnitRange()×BiInfUnitRange():\n   ⋮                       ⋮         ⋱  \n  30   25   20   15   10   5  0  -5  …  \n  24   20   16   12    8   4  0  -4     \n  18   15   12    9    6   3  0  -3     \n  12   10    8    6    4   2  0  -2     \n   6    5    4    3    2   1  0  -1     \n   0    0    0    0    0   0  0   0  …  \n  -6   -5   -4   -3   -2  -1  0   1     \n -12  -10   -8   -6   -4  -2  0   2     \n -18  -15  -12   -9   -6  -3  0   3     \n -24  -20  -16  -12   -8  -4  0   4     \n -30  -25  -20  -15  -10  -5  0   5  …  \n   ⋮                       ⋮         ⋱  " 
end