__precompile__()

module InfiniteArrays
   using Base, Compat, FillArrays

import Base: *, +, -, /, \, ==, isinf, isfinite, sign, angle, show, isless,
            fld, cld, div, min, max, minimum, maximum, mod,
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
            AbstractMatrix, AbstractArray, checkindex, unsafe_length, OneTo,
           to_shape, _sub2ind, print_matrix, print_matrix_row, print_matrix_vdots,
         checkindex, Slice, @_propagate_inbounds_meta, _in_range, _range, _rangestyle, Ordered,
         ArithmeticWraps, floatrange, reverse, unitrange_last,
         AbstractArray, AbstractVector
using Compat.LinearAlgebra
import Compat.LinearAlgebra: BlasInt, BlasFloat, norm, diag, diagm, ishermitian, issymmetric,
                             det, logdet, istriu, istril
import Compat: adjoint, axes


import Base: (:), _sub2ind_recurse
import LinearAlgebra: tr
const colon = (:)

export ∞, Hcat, Vcat, Zeros, Ones, Fill


include("Infinity.jl")
include("infrange.jl")
include("infdiagonal.jl")
include("inffill.jl")
include("infconcat.jl")
include("infarrayshow.jl")

end # module
