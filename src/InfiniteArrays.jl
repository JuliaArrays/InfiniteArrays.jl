__precompile__()

module InfiniteArrays
   using Base, Statistics, LinearAlgebra, FillArrays, LazyArrays, DSP

import Base: *, +, -, /, \, ==, isinf, isfinite, sign, angle, show, isless,
            fld, cld, div, min, max, minimum, maximum, mod,
            <, ≤, >, ≥, promote_rule, convert, copy,
            size, step, isempty, length, first, last,
            getindex, setindex!, intersect, @_inline_meta,
            sort, sort!, issorted, sortperm, sum, in, broadcast,
            eltype, parent, real, imag,
            conj, transpose,
            exp, log, sqrt, cos, sin, tan, csc, sec, cot,
            cosh, sinh, tanh, csch, sech, coth, acos, asin, atan, acsc, asec, acot,
            acosh, asinh, atanh, acsch, asech, acoth, (:),
            AbstractMatrix, AbstractArray, checkindex, unsafe_length, OneTo,
            to_shape, _sub2ind, print_matrix, print_matrix_row, print_matrix_vdots,
            checkindex, Slice, @propagate_inbounds, @_propagate_inbounds_meta,
         _in_range, _range, _rangestyle, Ordered,
         ArithmeticWraps, floatrange, reverse, unitrange_last,
         AbstractArray, AbstractVector, Array, Vector, Matrix,
         axes, (:), _sub2ind_recurse, broadcast, promote_eltypeof,
         diff, cumsum, show_delim_array, show_circular, Int,
         similar, _unsafe_getindex, string, zeros, fill, permutedims,
<<<<<<< Updated upstream
         cat_similar, vcat
=======
         cat_similar, vcat,
         reshape, ReshapedIndex, ind2sub_rs, _unsafe_getindex_rs
>>>>>>> Stashed changes

using Base.Broadcast
import Base.Broadcast: BroadcastStyle, AbstractArrayStyle, Broadcasted, broadcasted,
                        @nexprs, @ncall, combine_eltypes, DefaultArrayStyle, instantiate

import LinearAlgebra: BlasInt, BlasFloat, norm, diag, diagm, ishermitian, issymmetric,
                             det, logdet, istriu, istril, adjoint, tr, AbstractTriangular,
                             norm2, norm1, normp

import Statistics: mean, median

import FillArrays: AbstractFill, getindex_value
import LazyArrays: LazyArrayStyle, AbstractBandedLayout, MemoryLayout, LazyLayout,
                    ZerosLayout, @lazymul, AbstractArrayApplyStyle, CachedArray, CachedVector

import DSP: conv

export ∞, Hcat, Vcat, Zeros, Ones, Fill, Eye, BroadcastArray, cache





include("Infinity.jl")
include("infrange.jl")
include("infarrays.jl")

##
# Fill FillArrays
##

@lazymul Ones{<:Any,1,Tuple{OneToInf{Int}}}
@lazymul Fill{<:Any,1,Tuple{OneToInf{Int}}}
@lazymul Zeros{<:Any,1,Tuple{OneToInf{Int}}}

@lazymul Ones{<:Any,2,Tuple{OneToInf{Int},OneToInf{Int}}}
@lazymul Ones{<:Any,2,<:Tuple{OneToInf{Int},<:Any}}
@lazymul Ones{<:Any,2,<:Tuple{<:Any,OneToInf{Int}}}

@lazymul Fill{<:Any,2,Tuple{OneToInf{Int},OneToInf{Int}}}
@lazymul Fill{<:Any,2,<:Tuple{OneToInf{Int},<:Any}}
@lazymul Fill{<:Any,2,<:Tuple{<:Any,OneToInf{Int}}}

@lazymul Zeros{<:Any,2,Tuple{OneToInf{Int},OneToInf{Int}}}
@lazymul Zeros{<:Any,2,<:Tuple{OneToInf{Int},<:Any}}
@lazymul Zeros{<:Any,2,<:Tuple{<:Any,OneToInf{Int}}}

length(::Ones{<:Any,1,Tuple{OneToInf{Int}}}) = ∞
length(::Fill{<:Any,1,Tuple{OneToInf{Int}}}) = ∞
length(::Zeros{<:Any,1,Tuple{OneToInf{Int}}}) = ∞
length(::Ones{<:Any,2,Tuple{OneToInf{Int},OneToInf{Int}}}) = ∞
length(::Ones{<:Any,2,<:Tuple{OneToInf{Int},<:Any}}) = ∞
length(::Ones{<:Any,2,<:Tuple{<:Any,OneToInf{Int}}}) = ∞
length(::Fill{<:Any,2,Tuple{OneToInf{Int},OneToInf{Int}}}) = ∞
length(::Fill{<:Any,2,<:Tuple{OneToInf{Int},<:Any}}) = ∞
length(::Fill{<:Any,2,<:Tuple{<:Any,OneToInf{Int}}}) = ∞
length(::Zeros{<:Any,2,Tuple{OneToInf{Int},OneToInf{Int}}}) = ∞
length(::Zeros{<:Any,2,<:Tuple{OneToInf{Int},<:Any}}) = ∞
length(::Zeros{<:Any,2,<:Tuple{<:Any,OneToInf{Int}}}) = ∞

for op in (:norm2, :norm1)
   @eval $op(a::Zeros{T,N,NTuple{N,OneToInf{Int}}}) where {T,N} = norm(getindex_value(a))
end

normp(a::Zeros{T,N,NTuple{N,OneToInf{Int}}}, p) where {T,N} = norm(getindex_value(a))

for N=1:3
   for op in (:norm2, :norm1)
      @eval function $op(a::AbstractFill{T,$N,NTuple{$N,OneToInf{Int}}}) where {T,N}
         z = norm(getindex_value(a))
         iszero(z) && return z
         typeof(z)(Inf)
      end
   end
   @eval function normp(a::AbstractFill{T,$N,NTuple{$N,OneToInf{Int}}}, p) where {T,N }
      z = norm(getindex_value(a))
      iszero(z) && return z
      typeof(z)(Inf)
   end
end

for Typ in (:Number, :AbstractVector)
   @eval begin
      vcat(a::$Typ, b::AbstractFill{<:Any,1,<:Tuple{<:OneToInf}}) = Vcat(a, b)      
      vcat(a::$Typ, c::CachedVector{<:Any,<:Any,<:AbstractFill{<:Any,1,<:Tuple{<:OneToInf}}}) = 
         CachedArray(vcat(a, view(c.data,1:c.datasize[1])), c.array)
   end
end

# cat_similar(A, T, ::Tuple{Infinity}) = zeros(T, ∞)

##
# Temporary hacks for base support
##
OneTo(::Infinity) = OneToInf()
OneTo{T}(::Infinity) where T<:Integer = OneToInf{T}()
UnitRange(start::Integer, ::Infinity) = InfUnitRange(start)
UnitRange{T}(start::Integer, ::Infinity) where T<:Real = InfUnitRange{T}(start)
OneTo(a::OneToInf) = a
OneTo{T}(::OneToInf) where T<:Integer = OneToInf{T}()

Int(::Infinity) = ∞




end # module
