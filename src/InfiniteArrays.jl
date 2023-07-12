module InfiniteArrays
using ArrayLayouts: LayoutVecOrMat
using Base, Statistics, LinearAlgebra, FillArrays, Infinities, LazyArrays, ArrayLayouts

import Base: *, +, -, /, \, ==, isinf, isfinite, sign, signbit, angle, show, isless,
            fld, cld, div, min, max, minimum, maximum, mod,
            <, ≤, >, ≥, promote_rule, convert, unsafe_convert, copy,
            size, step, isempty, length, first, last, tail,
            getindex, setindex!, intersect, @_inline_meta,
            sort, sort!, issorted, sortperm, sum, in, broadcast,
            eltype, elsize, parent, parentindices, reinterpret,
            unaliascopy, dataids,
            real, imag, conj, transpose,
            exp, log, sqrt, cos, sin, tan, csc, sec, cot,
            cosh, sinh, tanh, csch, sech, coth, acos, asin, atan, acsc, asec, acot,
            acosh, asinh, atanh, acsch, asech, acoth, (:),
            AbstractMatrix, AbstractArray, checkindex, unsafe_length, unsafe_indices, OneTo,
            to_shape, _sub2ind, print_matrix, print_matrix_row, print_matrix_vdots,
            checkindex, Slice, IdentityUnitRange, @propagate_inbounds, @_propagate_inbounds_meta,
         	_in_range, _range, Ordered,
         	ArithmeticWraps, ArithmeticUnknown, floatrange, reverse, unitrange_last,
         	AbstractArray, AbstractVector, Array, Vector, Matrix,
         	axes, (:), _sub2ind_recurse, promote_eltypeof,
         	diff, cumsum, show_delim_array, show_circular, Int,
         	similar, _unsafe_getindex, string, zeros, ones, fill, permutedims,
         	cat_similar, vcat, hcat, one, zero,
		 	reshape, ReshapedIndex, ind2sub_rs, _unsafe_getindex_rs,
            searchsorted, searchsortedfirst, searchsortedlast, Ordering, lt, Fix2, findfirst,
            cat_indices, cat_size, cat_similar, __cat, _ind2sub_recurse, union, intersect, IEEEFloat

if VERSION < v"1.8-"
   import Base: _rangestyle
else
   import Base: range_start_step_length
end
if VERSION ≥ v"1.11.0-DEV.21"
   using LinearAlgebra: UpperOrLowerTriangular
else
   const UpperOrLowerTriangular{T,S} = Union{LinearAlgebra.UpperTriangular{T,S},
                                             LinearAlgebra.UnitUpperTriangular{T,S},
                                             LinearAlgebra.LowerTriangular{T,S},
                                             LinearAlgebra.UnitLowerTriangular{T,S}}
end


using Base.Broadcast
import Base.Broadcast: BroadcastStyle, AbstractArrayStyle, Broadcasted, broadcasted,
                        @nexprs, @ncall, combine_eltypes, DefaultArrayStyle, instantiate, axistype

import LinearAlgebra: BlasInt, BlasFloat, norm, diag, diagm, ishermitian, issymmetric,
                             det, logdet, istriu, istril, adjoint, tr, AbstractTriangular,
                             norm2, norm1, normp, AdjOrTrans, HermOrSym

import Statistics: mean, median

import FillArrays: AbstractFill, getindex_value, fill_reshape, RectDiagonal, Fill, Ones, Zeros, Eye
import LazyArrays: LazyArrayStyle, AbstractBandedLayout, MemoryLayout, LazyLayout, UnknownLayout,
                    ZerosLayout, AbstractCachedVector, CachedArray, CachedVector, ApplyLayout, LazyMatrix,
                    reshapedlayout, sub_materialize, sublayout, LayoutMatrix, LayoutVector, _padded_sub_materialize, PaddedLayout,
                    AbstractCachedMatrix, sub_paddeddata, InvColumnLayout

import ArrayLayouts: RangeCumsum, LayoutVecOrMat, LayoutVecOrMats
import Infinities: ∞, Infinity, InfiniteCardinal

export ∞, ℵ₀, Hcat, Vcat, Zeros, Ones, Fill, Eye, BroadcastArray, cache
import Base: unitrange, oneto



include("infrange.jl")
include("infarrays.jl")
include("reshapedarray.jl")

##
# Fill FillArrays
##

length(::Ones{<:Any,2,Tuple{OneToInf{Int},OneToInf{Int}}}) = ℵ₀
length(::Ones{<:Any,2,<:Tuple{OneToInf{Int},<:Any}}) = ℵ₀
length(::Ones{<:Any,2,<:Tuple{<:Any,OneToInf{Int}}}) = ℵ₀
length(::Fill{<:Any,2,Tuple{OneToInf{Int},OneToInf{Int}}}) = ℵ₀
length(::Fill{<:Any,2,<:Tuple{OneToInf{Int},<:Any}}) = ℵ₀
length(::Fill{<:Any,2,<:Tuple{<:Any,OneToInf{Int}}}) = ℵ₀
length(::Zeros{<:Any,2,Tuple{OneToInf{Int},OneToInf{Int}}}) = ℵ₀
length(::Zeros{<:Any,2,<:Tuple{OneToInf{Int},<:Any}}) = ℵ₀
length(::Zeros{<:Any,2,<:Tuple{<:Any,OneToInf{Int}}}) = ℵ₀

for op in (:norm2, :norm1)
   @eval $op(a::Zeros{T,N,NTuple{N,OneToInf{Int}}}) where {T,N} = norm(getindex_value(a))
end

normp(a::Zeros{T,N,NTuple{N,OneToInf{Int}}}, p) where {T,N} = norm(getindex_value(a))

for N=1:3
   for op in (:norm2, :norm1)
      @eval function $op(a::AbstractFill{T,$N,NTuple{$N,OneToInf{Int}}}) where T
         z = norm(getindex_value(a))
         iszero(z) && return z
         typeof(z)(Inf)
      end
   end
   @eval function normp(a::AbstractFill{T,$N,NTuple{$N,OneToInf{Int}}}, p) where T
      z = norm(getindex_value(a))
      iszero(z) && return z
      typeof(z)(Inf)
   end
end

for Typ in (:Number, :AbstractVector)
   @eval begin
      vcat(a::$Typ, b::AbstractFill{<:Any,1,<:Tuple{OneToInf}}) = Vcat(a, b)
      vcat(a::$Typ, c::CachedVector{<:Any,<:Any,<:AbstractFill{<:Any,1,<:Tuple{OneToInf}}}) =
         CachedArray(vcat(a, view(c.data,1:c.datasize[1])), c.array)
   end
end

vcat(a::AbstractVector, b::AbstractVector, c::AbstractFill{<:Any,1,<:Tuple{OneToInf}}) = Vcat(vcat(a,b), c)
vcat(a::AbstractVector, b::AbstractVector, c::AbstractVector, d::AbstractFill{<:Any,1,<:Tuple{OneToInf}}) = Vcat(vcat(a,b,c), d)

vcat(a::AbstractMatrix, b::AbstractFill{<:Any,2,<:Tuple{OneToInf,OneTo}}) = Vcat(a, b)

cat_similar(A, ::Type{T}, shape::Tuple{PosInfinity}) where T = zeros(T,∞)
cat_similar(A::AbstractArray, ::Type{T}, shape::Tuple{PosInfinity}) where T = zeros(T,∞)
cat_similar(A::Array, ::Type{T}, shape::Tuple{PosInfinity}) where T = zeros(T,∞)
function Base.__cat(A, shape::NTuple{N,PosInfinity}, catdims, X...) where N
   offsets = zeros(Union{Int,InfiniteCardinal{0}}, N)
   inds = Vector{Union{UnitRange{Int},InfUnitRange{Int}}}(undef, N)
   concat = copyto!(zeros(Bool, N), catdims)
   for x in X
       for i = 1:N
           if concat[i]
               inds[i] = offsets[i] .+ cat_indices(x, i)
               offsets[i] += cat_size(x, i)
           else
               inds[i] = 1:shape[i]
           end
       end
       I::NTuple{N, Union{InfUnitRange{Int},UnitRange{Int}}} = (inds...,)
       if x isa AbstractArray
           copyto!(view(A, I...), x)
       else
           fill!(view(A, I...), x)
       end
   end
   return A
end

reshape(parent::AbstractArray, shp::Tuple{OneToInf, Vararg{Union{Integer,OneTo,OneToInf}}}) =
   reshape(parent, to_shape(shp))
reshape(parent::AbstractArray, shp::Tuple{Union{Integer,OneTo}, OneToInf, Vararg{Union{Integer,OneTo,OneToInf}}}) =
   reshape(parent, to_shape(shp))


# cat_similar(A, T, ::Tuple{PosInfinity}) = zeros(T, ∞)

axistype(::OneTo{T}, ::OneToInf{V}) where {T,V} = OneToInf{promote_type(T,V)}()
axistype(::OneToInf{V}, ::OneTo{T}) where {T,V} = OneToInf{promote_type(T,V)}()

# sort.jl
# returns the range of indices of v equal to x
# if v does not contain x, returns a 0-length range
# indicating the insertion point of x
function searchsorted(v::AbstractVector, x, ilo::Int, ::PosInfinity, o::Ordering)
    lo = ilo-1
    hi = ℵ₀
    @inbounds while lo < hi-1
        m = isinf(hi) ? lo + 1000 : (lo+hi)>>>1
        if lt(o, v[m], x)
            lo = m
        elseif lt(o, x, v[m])
            hi = m
        else
            a = searchsortedfirst(v, x, max(lo,ilo), m, o)
            b = searchsortedlast(v, x, m, hi, o)
            return a : b
        end
    end
    return (lo + 1) : (hi - 1)
end


# index of the first value of vector a that is greater than or equal to x;
# returns length(v)+1 if x is greater than all values in v.
function searchsortedfirst(v::AbstractVector, x, lo::Int, hi::PosInfinity, o::Ordering)
   u = 1
   lo = lo - u
   hi = ℵ₀
   @inbounds while lo < hi - u
      m = isinf(hi) ? lo + 1000 : (lo+hi)>>>1
      if lt(o, v[m], x)
         lo = m
      else
         hi = m
      end
   end
   return hi
end

# index of the last value of vector a that is less than or equal to x;
# returns 0 if x is less than all values of v.
function searchsortedlast(v::AbstractVector, x, lo::Int, hi::PosInfinity, o::Ordering)
   u = 1
   lo = lo - u
   hi = ℵ₀
   @inbounds while lo < hi - u
       m = isinf(hi) ? lo + 1000 : (lo+hi)>>>1
       if lt(o, x, v[m])
           hi = m
       else
           lo = m
       end
   end
   return lo
end

# special case for Vcat
@inline function LazyArrays.searchsortedlast_recursive(::PosInfinity, x, a, args...)
    n = sum(map(length,args))
    r = searchsortedlast(a, x)
    r > 0 && return n + r
    return LazyArrays.searchsortedlast_recursive(n, x, args...)
end




end # module
