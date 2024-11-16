module InfiniteArrays
using LinearAlgebra, FillArrays, Infinities, LazyArrays, ArrayLayouts

import Base: *, +, -, /, <, ==, >, \, ≤, ≥, (:), @propagate_inbounds,
             AbstractArray, AbstractMatrix, AbstractVector, ArithmeticUnknown, ArithmeticWraps, Array, Fix2, IEEEFloat,
             IdentityUnitRange, Int, Matrix, OneTo, Ordered, Ordering, ReshapedIndex, Slice, Vector, __cat, _in_range,
             _ind2sub_recurse, _range, _sub2ind, _sub2ind_recurse, _unsafe_getindex, _unsafe_getindex_rs,
             angle, axes, broadcast, cat_indices,
             cat_similar, cat_size, checkindex, collect, convert, copy,
             cumsum, dataids, diff, div, eltype, fill, findfirst, first, floatrange, getindex, hcat,
             in, ind2sub_rs, intersect, inv, isempty, isinf, issorted, last, length, lt, max,
             maximum, minimum, mod, one, ones, parent, parentindices, permutedims, print_matrix, print_matrix_row,
             print_matrix_vdots, promote_rule, reinterpret, reshape, reverse, searchsorted,
             searchsortedfirst, searchsortedlast, setindex!, show, show_circular, show_delim_array, sign,
             signbit, similar, size, sort, sort!, step, sum, tail,
             to_shape, transpose, unaliascopy, union, unitrange_last, unsafe_convert, unsafe_indices, unsafe_length,
             vcat, zeros

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
import ArrayLayouts: AbstractBandedLayout, LayoutMatrix, LayoutVecOrMat, LayoutVecOrMats, LayoutVector, MemoryLayout,
                     RangeCumsum, UnknownLayout, reshapedlayout, sub_materialize, sublayout

import Base.Broadcast: BroadcastStyle, Broadcasted, DefaultArrayStyle, axistype, broadcasted

import FillArrays: AbstractFill, getindex_value, fill_reshape, RectDiagonal, Fill, Ones, Zeros, Eye, elconvert

import Infinities: InfiniteCardinal, Infinity, ∞

import LazyArrays: AbstractCachedVector, ApplyLayout, CachedArray, CachedVector, InvColumnLayout,
                   LazyArrayStyle, LazyLayout, LazyMatrix, PaddedColumns, _padded_sub_materialize, sub_paddeddata

import LinearAlgebra: AdjOrTrans, HermOrSym, diag, norm, norm1, norm2, normp

import LazyArrays: AbstractPaddedLayout

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

for Typ in (:Number, :(AbstractVector{<:Number}))
   @eval begin
      vcat(a::$Typ, b::AbstractFill{<:Number,1,<:Tuple{OneToInf}}) = Vcat(a, b)
      vcat(a::$Typ, c::CachedVector{<:Number,<:Any,<:AbstractFill{<:Any,1,<:Tuple{OneToInf}}}) =
         CachedArray(vcat(a, view(c.data,1:c.datasize[1])), c.array)
   end
end

vcat(a::AbstractVector{<:Number}, b::AbstractVector{<:Number}, c::AbstractFill{<:Number,1,<:Tuple{OneToInf}}) = Vcat(vcat(a,b), c)
vcat(a::AbstractVector{<:Number}, b::AbstractVector{<:Number}, c::AbstractVector{<:Number}, d::AbstractFill{<:Number,1,<:Tuple{OneToInf}}) = Vcat(vcat(a,b,c), d)

vcat(a::AbstractMatrix{<:Number}, b::AbstractFill{<:Number,2,<:Tuple{OneToInf,OneTo}}) = Vcat(a, b)

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
function searchsorted(v::AbstractVector, x, ilo::Integer, ::PosInfinity, o::Ordering)
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
function searchsortedfirst(v::AbstractVector, x, lo::Integer, hi::PosInfinity, o::Ordering)
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
function searchsortedlast(v::AbstractVector, x, lo::Integer, hi::PosInfinity, o::Ordering)
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


collect(G::Base.Generator{<:InfRanges}) = BroadcastArray(G.f, G.iter)

if !isdefined(Base, :get_extension)
    include("../ext/InfiniteArraysStatisticsExt.jl")
    include("../ext/InfiniteArraysDSPExt.jl")
end


end # module
