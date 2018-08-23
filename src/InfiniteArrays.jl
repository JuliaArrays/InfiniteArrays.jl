__precompile__()

module InfiniteArrays
   using Base, Statistics, LinearAlgebra, FillArrays, LazyArrays

import Base: *, +, -, /, \, ==, isinf, isfinite, sign, angle, show, isless,
            fld, cld, div, min, max, minimum, maximum, mod,
            <, ≤, >, ≥, promote_rule, convert, copy,
            size, step, isempty, length, first, last, 
            getindex, setindex!, OneTo, intersect, @_inline_meta,
            sort, sort!, issorted, sortperm, sum, in, broadcast,
            eltype, parent, real, imag,
            conj, transpose,
            exp, log, sqrt,
                      cos, sin, tan, csc, sec, cot,
                      cosh, sinh, tanh, csch, sech, coth,
                      acos, asin, atan, acsc, asec, acot,
                      acosh, asinh, atanh, acsch, asech, acoth, (:),
            AbstractMatrix, AbstractArray, checkindex, unsafe_length, OneTo,
           to_shape, _sub2ind, print_matrix, print_matrix_row, print_matrix_vdots,
         checkindex, Slice, @propagate_inbounds, @_propagate_inbounds_meta,
         _in_range, _range, _rangestyle, Ordered,
         ArithmeticWraps, floatrange, reverse, unitrange_last,
         AbstractArray, AbstractVector, axes, (:), _sub2ind_recurse, broadcast, promote_eltypeof

using Base.Broadcast
import Base.Broadcast: BroadcastStyle, AbstractArrayStyle, Broadcasted, broadcasted,
                        @nexprs, @ncall, combine_eltypes, DefaultArrayStyle, instantiate

import LinearAlgebra: BlasInt, BlasFloat, norm, diag, diagm, ishermitian, issymmetric,
                             det, logdet, istriu, istril, adjoint, tr

import Statistics: mean, median

import LazyArrays: LazyArrayStyle

export ∞, Hcat, Vcat, Zeros, Ones, Fill, Eye, BroadcastArray, cache


include("Infinity.jl")
include("infrange.jl")
include("infdiagonal.jl")
include("inffill.jl")
include("infarrayshow.jl")

end # module
