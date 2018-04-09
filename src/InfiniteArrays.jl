__precompile__()

module InfiniteArrays
   using Base, Compat, FillArrays

import Base: *, +, -, /, \, ==, isinf, isfinite, sign, angle, show, isless,
            fld, cld, div, min, max, minimum, maximum,
            <, ≤, >, ≥, promote_rule, convert,
            size, step, isempty, length, first, last, start, next, done,
            getindex, setindex!, OneTo, intersect, @_inline_meta,
            sort, sort!, issorted, sortperm, sum, mean, median, in, broadcast,
            eltype, parent, real, imag,
            conj, transpose,
            exp, log, sqrt,
                      cos, sin, tan, csc, sec, cot,
                      cosh, sinh, tanh, csch, sech, coth,
                      acos, asin, atan, acsc, asec, acot,
                      acosh, asinh, atanh, acsch, asech, acoth, (:),
            AbstractMatrix, AbstractArray, axes, inds2string, checkindex, unsafe_length, OneTo,
           to_shape
using Compat.LinearAlgebra
import Compat.LinearAlgebra: BlasInt, BlasFloat, norm, diag, diagm, ishermitian, issymmetric,
                             det, logdet, istriu, istril, tr
import Compat: adjoint, axes

if VERSION ≥ v"0.7-"
   import Base: (:)
   const colon = (:)
else
   import Base: colon
end

export ∞

abstract type InfArray{T,N} <: AbstractArray{T,N} end
const InfVector{T} = InfArray{T,1}
const InfMatrix{T} = InfArray{T,2}

function axes(A::InfArray)
    @_inline_meta
    map(OneTo, size(A))
end


include("Infinity.jl")

include("infrange.jl")
include("infdiagonal.jl")
include("inffill.jl")

end # module
