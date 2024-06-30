"""
    InfRandVector([rng=default_rng()], [dist=Float64])

Represents a random infinite sequence. The random number generator can be specified 
by the first argument `rng`, which defaults to `Random.default_rng()`, and the distribution 
to generate from can be specified using the `dist` argument, which defaults to `Float64`.

```julia-repl
julia> using InfiniteArrays

julia> seq = InfiniteArrays.InfRandVector();

julia> seq[1]
0.6947847847152157

julia> seq[1:5]
5-element Vector{Float64}:
 0.6947847847152157
 0.49859004150722164
 0.31559745572937126
 0.5338596092137163
 0.14792462133894646

julia> using Distributions, Random

julia> seq = InfiniteArrays.InfRandVector(MersenneTwister(123), Normal(0.3, 1.7))
ℵ₀-element InfRandVector{Float64, Normal{Float64}, MersenneTwister} with indices 1:∞:
  2.3234553976766703
  3.7819055032417075
  2.242506534874238
  1.0810065546920364
 -0.3743544348018791
 -0.8300113268258689
  1.967645305489507
  ⋮
```
"""
mutable struct InfRandVector{T,D,RNG} <: AbstractCachedVector{T}
    const rng::RNG
    const dist::D
    const data::Vector{T}
    datasize::Int
end
function InfRandVector(rng=default_rng(), dist=Float64)
    T = typeof(rand(copy(rng), dist))
    _rng = copy(rng)
    return InfRandVector{T,typeof(dist),typeof(_rng)}(_rng, dist, T[], 0)
end
size(::InfRandVector) = (ℵ₀,)
axes(::InfRandVector) = (1:ℵ₀,)
length(::InfRandVector) = ℵ₀
function resizedata!(seq::InfRandVector, inds)
    newlen = maximum(inds)
    curlen = length(seq.data)
    newlen > curlen || return seq
    resize!(seq.data, newlen)
    # rand!(seq.rng, view(seq.data, curlen+1:newlen), seq.dist)
    # ^ rand() is not actually sequential.. rand(Random.seed!(123), 1000) ≠ (rng = Random.seed!(123); [rand(rng) for _ in 1:1000])
    for i in (curlen+1):newlen
        seq.data[i] = rand(seq.rng, seq.dist)
    end
    seq.datasize = newlen
    return seq
end

"""
    InfRandMatrix([rng=default_rng()], n; dist=Float64])

Represents a random infinite matrix with `n` rows. The random number generator 
can be specified by the first argument `rng`, which defaults to `Random.default_rng()`, and the number of 
rows is given by the `n` argument. The `dist` keyword argument (default `Float64`)
can be used to specify the distribution to sample from. 
"""
struct InfRandMatrix{T, S <: InfRandVector{T}} <: LazyMatrix{T} 
    seq::S 
    n::Int
end 
InfRandMatrix(n::Int; dist=Float64) = InfRandMatrix(default_rng(), n; dist)
InfRandMatrix(rng, n::Int; dist=Float64) = InfRandMatrix(InfRandVector(rng, dist), n)
function Base.getindex(A::InfRandMatrix, i::Int, j::Int)
    ((i < 1) || (i > A.n) || (j < 0)) && throw(BoundsError(A, (i, j)))
    lin = (j - 1) * A.n + i 
    return A.seq[lin]
end
size(A::InfRandMatrix) = (A.n, ∞)

