# This file is based on part of Julia v0.7. License is MIT: https://julialang.org/license

## Diagonal matrices

struct InfDiagonal{T,V<:InfVector} <: InfMatrix{T}
    diag::V
    k::Int
end
"""
    InfDiagonal(A::InfMatrix)

Construct a matrix from the diagonal of `A`.
"""
InfDiagonal(A::InfMatrix) = InfDiagonal(diag(A))

"""
    diagm(V::InfVector)

Construct a matrix with `V` as its diagonal.

# Examples
```jldoctest
julia> V = [1, 2]
2-element Array{Int64,1}:
 1
 2

julia> diagm(V)
2×2 Diagonal{Int64,Array{Int64,1}}:
 1  ⋅
 ⋅  2
```
"""
InfDiagonal(V::InfVector{T}, k::Integer=0) where {T} = InfDiagonal{T,typeof(V)}(V,k)
InfDiagonal{T}(V::InfVector{T}, k::Integer=0) where {T} = InfDiagonal{T,typeof(V)}(V,k)
InfDiagonal{T}(V::InfVector, k::Integer=0) where {T} = InfDiagonal{T}(convert(InfVector{T}, V),k)

diagm(V::InfVector{T}, k::Integer=0) where T = InfDiagonal{T}(V,k)

InfDiagonal{T}(D::InfDiagonal{T}) where {T} = D
InfDiagonal{T}(D::InfDiagonal) where {T} = InfDiagonal{T}(convert(InfVector{T}, D.diag))
AbstractMatrix{T}(D::InfDiagonal) where {T} = InfDiagonal{T}(D)

size(D::InfDiagonal) = (length(D.diag),length(D.diag))

function size(D::InfDiagonal,d::Integer)
    if d<1 || d > 2
        throw(ArgumentError("dimension must be ≥ 1, got $d"))
    end
    return length(D.diag)
end

@inline function getindex(D::InfDiagonal, i::Int, j::Int)
    @boundscheck checkbounds(D, i, j)
    if i - j == D.k
        @inbounds r = D.diag[min(i,j)]
    else
        r = diagzero(D, i, j)
    end
    r
end
diagzero(::InfDiagonal{T},i,j) where {T} = zero(T)
diagzero(D::InfDiagonal{Matrix{T}},i,j) where {T} = zeros(T, size(D.diag[i], 1), size(D.diag[j], 2))

function setindex!(D::InfDiagonal, v, i::Int, j::Int)
    @boundscheck checkbounds(D, i, j)
    if i == j
        @inbounds D.diag[i] = v
    elseif !iszero(v)
        throw(ArgumentError("cannot set off-diagonal entry ($i, $j) to a nonzero value ($v)"))
    end
    return v
end


## structured matrix methods ##
parent(D::InfDiagonal) = D.diag

ishermitian(D::InfDiagonal{<:Real}) = true
ishermitian(D::InfDiagonal{<:Number}) = isreal(D.diag)
issymmetric(D::InfDiagonal{<:Number}) = true


factorize(D::InfDiagonal) = D

broadcast(::typeof(abs), D::InfDiagonal) = InfDiagonal(abs.(D.diag))
real(D::InfDiagonal) = InfDiagonal(real(D.diag))
imag(D::InfDiagonal) = InfDiagonal(imag(D.diag))

istriu(D::InfDiagonal) = D.k ≥ 0
istril(D::InfDiagonal) = D.k ≤ 0

(==)(Da::InfDiagonal, Db::InfDiagonal) = Da.k == Db.k && Da.diag == Db.diag
(-)(A::InfDiagonal) = InfDiagonal(-A.diag)
# (+)(Da::InfDiagonal, Db::InfDiagonal) = InfDiagonal(Da.diag + Db.diag)
# (-)(Da::InfDiagonal, Db::InfDiagonal) = InfDiagonal(Da.diag - Db.diag)

(*)(x::Number, D::InfDiagonal) = InfDiagonal(x * D.diag)
(*)(D::InfDiagonal, x::Number) = InfDiagonal(D.diag * x)
(/)(D::InfDiagonal, x::Number) = InfDiagonal(D.diag / x)
# (*)(Da::InfDiagonal, Db::InfDiagonal) = InfDiagonal(Da.diag .* Db.diag)
function (*)(D::InfDiagonal, V::InfVector)
    if D.k == 0
        D.diag .* V
    elseif D.k > 0
        D.diag .* V[D.k:end]
    else # D.k < 0
        [zeros(-D.k); D.diag .* V]
    end
end

conj(D::InfDiagonal) = InfDiagonal(conj(D.diag))
transpose(D::InfDiagonal{<:Number}) = D
transpose(D::InfDiagonal) = InfDiagonal(transpose.(D.diag))
adjoint(D::InfDiagonal{<:Number}) = conj(D)
adjoint(D::InfDiagonal) = InfDiagonal(adjoint.(D.diag))

function diag(D::InfDiagonal, k::Integer=0)
    # every branch call similar(..., ::Int) to make sure the
    # same vector type is returned independent of k
    if k == D.k
        return D.diag
    else
        return zeros(∞)
    end
end
tr(D::InfDiagonal{T}) where T = D.k == 0 ? zero(T) : T(sum(D.diag))
det(D::InfDiagonal{T}) where T = D.k == 0 ? zero(T) : T(prod(D.diag))
logdet(D::InfDiagonal{T})  where T<:Real = D.k == 0 ? log(zero(T)) : T(sum(log, D.diag))
function logdet(D::InfDiagonal{<:Complex}) # make sure branch cut is correct
    z = sum(log, D.diag)
    complex(real(z), rem2pi(imag(z), RoundNearest))
end

# Matrix functions
for f in (:exp, :log, :sqrt,
          :cos, :sin, :tan, :csc, :sec, :cot,
          :cosh, :sinh, :tanh, :csch, :sech, :coth,
          :acos, :asin, :atan, :acsc, :asec, :acot,
          :acosh, :asinh, :atanh, :acsch, :asech, :acoth)
    @eval $f(D::InfDiagonal) = InfDiagonal($f.(D.diag))
end
