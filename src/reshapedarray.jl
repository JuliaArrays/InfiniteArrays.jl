# This file is based on a part of Julia. License is MIT: https://julialang.org/license

using  Base.MultiplicativeInverses: SignedMultiplicativeInverse

struct ReshapedArray{T,N,P<:AbstractArray,DIMS<:Tuple,MI<:Tuple{Vararg{SignedMultiplicativeInverse{Int}}}} <: AbstractArray{T,N}
    parent::P
    dims::DIMS
    mi::MI
end
ReshapedArray(parent::AbstractArray{T}, dims::NTuple{N,Integer}, mi) where {T,N} = ReshapedArray{T,N,typeof(parent),typeof(dims),typeof(mi)}(parent, dims, mi)
ReshapedArray(parent::AbstractArray{T}, dims::NTuple{N,Integer}) where {T,N} = ReshapedArray(parent, dims, ())

size(A::ReshapedArray) = A.dims
similar(A::ReshapedArray, eltype::Type, dims::Dims) = similar(parent(A), eltype, dims)
parent(A::ReshapedArray) = A.parent
parentindices(A::ReshapedArray) = map(oneto, size(parent(A)))
# reinterpret(::Type{T}, A::ReshapedArray, dims::Dims) where {T} = reinterpret(T, parent(A), dims)

unaliascopy(A::ReshapedArray) = typeof(A)(unaliascopy(A.parent), A.dims, A.mi)
dataids(A::ReshapedArray) = dataids(A.parent)

@inline function getindex(A::ReshapedArray{T,N}, indices::Vararg{Integer,N}) where {T,N}
    @boundscheck checkbounds(A, indices...)
    _unsafe_getindex(A, indices...)
end
@inline function getindex(A::ReshapedArray, index::ReshapedIndex)
    @boundscheck checkbounds(parent(A), index.parentindex)
    @inbounds ret = parent(A)[index.parentindex]
    ret
end

@inline function _unsafe_getindex(A::ReshapedArray{T,N}, indices::Vararg{Integer,N}) where {T,N}
    i = Base._sub2ind(size(A), indices...)
    I = ind2sub_rs(axes(A.parent), A.mi, i)
    _unsafe_getindex_rs(parent(A), I)
end

@inline function setindex!(A::ReshapedArray{T,N}, val, indices::Vararg{Integer,N}) where {T,N}
    @boundscheck checkbounds(A, indices...)
    _unsafe_setindex!(A, val, indices...)
end
@inline function setindex!(A::ReshapedArray, val, index::ReshapedIndex)
    @boundscheck checkbounds(parent(A), index.parentindex)
    @inbounds parent(A)[index.parentindex] = val
    val
end

@inline function _unsafe_setindex!(A::ReshapedArray{T,N}, val, indices::Vararg{Integer,N}) where {T,N}
    @inbounds parent(A)[ind2sub_rs(axes(A.parent), A.mi, Base._sub2ind(size(A), indices...))...] = val
    val
end

unsafe_convert(::Type{Ptr{T}}, a::ReshapedArray{T}) where {T} = unsafe_convert(Ptr{T}, parent(a))

reshape(A::AbstractArray, a::PosInfinity, b::Integer...) = ReshapedArray(A, tuple(ℵ₀,b...))
reshape(A::AbstractArray, a::Integer, b::PosInfinity, c::Integer...) = ReshapedArray(A, tuple(a,ℵ₀,c...))
reshape(A::AbstractArray, dims::Tuple{PosInfinity,Vararg{Integer}}) = ReshapedArray(A, Base.to_shape(dims))
reshape(A::AbstractArray, dims::Tuple{Integer,PosInfinity,Vararg{Integer}}) = ReshapedArray(A, Base.to_shape(dims))
reshape(A::AbstractFill, a::PosInfinity, b::Integer...) = fill_reshape(A, ℵ₀, b...)
reshape(A::AbstractFill, a::Integer, b::PosInfinity, c::Integer...) = fill_reshape(A, a, ℵ₀, c...)
reshape(A::AbstractFill, dims::Tuple{PosInfinity,Vararg{Integer}}) = fill_reshape(A, Base.to_shape(dims)...)
reshape(A::AbstractFill, dims::Tuple{Integer,PosInfinity,Vararg{Integer}}) = fill_reshape(A, Base.to_shape(dims)...)


BroadcastStyle(::Type{<:ReshapedArray{T,N,<:Any,NTuple{N,InfiniteCardinal{0}}}}) where {T,N} = LazyArrayStyle{N}()
BroadcastStyle(::Type{<:ReshapedArray{T,2,<:Any,<:Tuple{Any,InfiniteCardinal{0}}}}) where {T} = LazyArrayStyle{2}()
BroadcastStyle(::Type{<:ReshapedArray{T,2,<:Any,<:Tuple{InfiniteCardinal{0},Any}}}) where {T} = LazyArrayStyle{2}()


MemoryLayout(::Type{<:ReshapedArray{T,N,A,DIMS}}) where {T,N,A,DIMS} = reshapedlayout(MemoryLayout(A), DIMS)


###
# permutedims for reshaped unrolls
###

permutedims(R::ReshapedArray{<:Any,2,<:AbstractVector}) = parent(R)



###
# support Reshaping infinite-vector to matrix
function Base._sub2ind_recurse(inds, L::InfiniteCardinal{0}, ind, i::Integer, I::Integer...)
    r1 = inds[1]
    @assert iszero(Base.offsetin(i, r1))
    _sub2ind_recurse(tail(inds), Base.nextL(L, r1), ind, I...)
end
