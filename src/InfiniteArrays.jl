__precompile__()

module InfiniteArrays
   using Base, Compat

import Base: *, +, -, /, \, ==, colon, isinf, isfinite, sign, angle, show, isless,
            fld, cld, div, min, max, minimum, maximum,
            <, ≤, >, ≥, promote_rule, convert,
            size, step, isempty, length, first, last, start, next, done,
            getindex, OneTo, intersect, @_inline_meta,
            sort, sort!, issorted, sortperm, sum, mean, median, in, broadcast,
            eltype
import Compat.LinearAlgebra: BlasInt, BlasFloat, norm, (:)


export ∞

abstract type InfiniteArray{T,N} end
const InfiniteVector{T} = InfiniteArray{T,1}
const InfiniteMatrix{T} = InfiniteArray{T,2}

eltype(::Type{IA}) where IA<:InfiniteArray{T} where T = T
eltype(::IA) where IA<:InfiniteArray = eltype(IA)




include("Infinity.jl")
include("infiniterange.jl")

end # module
