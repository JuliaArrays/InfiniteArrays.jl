__precompile__()

module InfiniteArrays
   using Base, Compat

import Base: *, +, -, /, \, ==, colon, isinf, isfinite, sign, angle, show, isless,
            fld, cld, div, min, max, minimum, maximum,
            <, ≤, >, ≥, promote_rule, convert,
            size, step, isempty, length, first, last, start, next, done,
            getindex, setindex!, OneTo, intersect, @_inline_meta,
            sort, sort!, issorted, sortperm, sum, mean, median, in, broadcast,
            eltype, parent, ishermitian, issymmetric, real, imag, istriu, istril,
            conj, transpose, det, logdet,
            exp, log, sqrt,
                      cos, sin, tan, csc, sec, cot,
                      cosh, sinh, tanh, csch, sech, coth,
                      acos, asin, atan, acsc, asec, acot,
                      acosh, asinh, atanh, acsch, asech, acoth, (:),
            AbstractMatrix, AbstractArray, indices, inds2string, diag, diagm
import Compat.LinearAlgebra: BlasInt, BlasFloat, norm
import Compat: adjoint

if VERSION ≥ v"0.7-"
   import Base: tr
end

export ∞

abstract type InfArray{T,N}  end  # <: AbstractArray{T,N}
const InfVector{T} = InfArray{T,1}
const InfMatrix{T} = InfArray{T,2}

eltype(::Type{IA}) where IA<:InfArray{T} where T = T
eltype(::IA) where IA<:InfArray = eltype(IA)


oneto(n::Integer) = Base.OneTo(n)

function indices(A::InfArray)
    @_inline_meta
    map(oneto, size(A))
end


include("Infinity.jl")
include("infrange.jl")
include("infdiagonal.jl")

end # module
