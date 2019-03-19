# InfiniteArrays.jl

[![Travis](https://travis-ci.org/JuliaArrays/InfiniteArrays.jl.svg?branch=master)](https://travis-ci.org/JuliaArrays/InfiniteArrays.jl)
[![codecov](https://codecov.io/gh/JuliaArrays/InfiniteArrays.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaArrays/InfiniteArrays.jl)


A Julia package for representing infinite-dimensional arrays, designed to work
with other array types. Infinite arrays are by necessity lazy, and so this
package is closely linked to [LazyArrays.jl](https://github.com/JuliaArrays/LazyArrays.jl).

The package currently supports construction, but future versions will include
support for linear algebra operations like matrix multiplication, LU decompositions,
and QR decompositions.


Here are some examples of current functionality:

```julia
julia> using InfiniteArrays, LinearAlgebra

julia> x = Ones(∞) # infinite vector of all ones
Ones{Float64,1,Tuple{InfiniteArrays.OneToInf{Int64}}} with indices OneToInf():
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 ⋮

julia> cumsum(x) # infinite sum
1.0:1.0:+∞

julia> exp.(-(1:∞)) .+ 2 # broadcasting supported
BroadcastArray{Float64,1,Base.Broadcast.Broadcasted{LazyArrays.LazyArrayStyle{1},Tuple{InfiniteArrays.OneToInf{Int64}},typeof(+),Tuple{BroadcastArray{Float64,1,Base.Broadcast.Broadcasted{LazyArrays.LazyArrayStyle{1},Tuple{InfiniteArrays.OneToInf{Int64}},typeof(exp),Tuple{InfiniteArrays.InfStepRange{Int64,Int64}}}},Int64}}} with indices OneToInf():
 2.3678794411714423
 2.135335283236613
 2.049787068367864
 2.018315638888734
 2.0067379469990856
 2.0024787521766663
 2.0009118819655547
 2.0003354626279024
 2.0001234098040865
 2.0000453999297623
 2.0000167017007904
 2.000006144212353
 2.000002260329407
 2.000000831528719
 2.0000003059023204
 2.0000001125351745
 2.0000000413993773
 2.00000001522998
 2.0000000056027964
 2.0000000020611535
 2.000000000758256
 2.0000000002789466
 2.000000000102619
 2.000000000037751
 2.000000000013888
 2.0000000000051092
 ⋮

julia> Diagonal(1:∞) # combines well with Base arrays
Diagonal{Int64,InfiniteArrays.InfUnitRange{Int64}} with indices OneToInf()×OneToInf():
 1  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅   ⋅   ⋅   ⋅   ⋅   ⋅   ⋅  …
 ⋅  2  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅   ⋅   ⋅   ⋅   ⋅   ⋅   ⋅
 ⋅  ⋅  3  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅   ⋅   ⋅   ⋅   ⋅   ⋅   ⋅
 ⋅  ⋅  ⋅  4  ⋅  ⋅  ⋅  ⋅  ⋅   ⋅   ⋅   ⋅   ⋅   ⋅   ⋅
 ⋅  ⋅  ⋅  ⋅  5  ⋅  ⋅  ⋅  ⋅   ⋅   ⋅   ⋅   ⋅   ⋅   ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  6  ⋅  ⋅  ⋅   ⋅   ⋅   ⋅   ⋅   ⋅   ⋅  …
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  7  ⋅  ⋅   ⋅   ⋅   ⋅   ⋅   ⋅   ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  8  ⋅   ⋅   ⋅   ⋅   ⋅   ⋅   ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  9   ⋅   ⋅   ⋅   ⋅   ⋅   ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  10   ⋅   ⋅   ⋅   ⋅   ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅   ⋅  11   ⋅   ⋅   ⋅   ⋅  …
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅   ⋅   ⋅  12   ⋅   ⋅   ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅   ⋅   ⋅   ⋅  13   ⋅   ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅   ⋅   ⋅   ⋅   ⋅  14   ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅   ⋅   ⋅   ⋅   ⋅   ⋅  15
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅   ⋅   ⋅   ⋅   ⋅   ⋅   ⋅  …
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅   ⋅   ⋅   ⋅   ⋅   ⋅   ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅   ⋅   ⋅   ⋅   ⋅   ⋅   ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅   ⋅   ⋅   ⋅   ⋅   ⋅   ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅   ⋅   ⋅   ⋅   ⋅   ⋅   ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅   ⋅   ⋅   ⋅   ⋅   ⋅   ⋅  …
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅   ⋅   ⋅   ⋅   ⋅   ⋅   ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅   ⋅   ⋅   ⋅   ⋅   ⋅   ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅   ⋅   ⋅   ⋅   ⋅   ⋅   ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅   ⋅   ⋅   ⋅   ⋅   ⋅   ⋅
 ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅   ⋅   ⋅   ⋅   ⋅   ⋅   ⋅  …
 ⋮              ⋮                ⋮                  ⋱

 julia> C = cache(Ones(∞)); C[1] = 3; C # Use cache to make mutable
LazyArrays.CachedArray{Float64,1,Array{Float64,1},Ones{Float64,1,Tuple{InfiniteArrays.Infinity}}} with indices OneToInf():
 3.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 ⋮
 ```
