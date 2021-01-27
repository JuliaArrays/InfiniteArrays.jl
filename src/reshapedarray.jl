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
parentindices(A::ReshapedArray) = map(OneTo, size(parent(A)))
reinterpret(::Type{T}, A::ReshapedArray, dims::Dims) where {T} = reinterpret(T, parent(A), dims)
elsize(::Type{<:ReshapedArray{<:Any,<:Any,P}}) where {P} = elsize(P)

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

reshape(A::AbstractArray, a::Infinity, b::Integer...) = ReshapedArray(A, tuple(a,b...))
reshape(A::AbstractArray, a::Integer, b::Infinity, c::Integer...) = ReshapedArray(A, tuple(a,b,c...))
reshape(A::AbstractArray, dims::Tuple{Infinity,Vararg{Integer}}) = ReshapedArray(A, dims)
reshape(A::AbstractArray, dims::Tuple{Integer,Infinity,Vararg{Integer}}) = ReshapedArray(A, dims)
reshape(A::AbstractFill, a::Infinity, b::Integer...) = fill_reshape(A, a, b...)
reshape(A::AbstractFill, a::Integer, b::Infinity, c::Integer...) = fill_reshape(A, a, b, c...)
reshape(A::AbstractFill, dims::Tuple{Infinity,Vararg{Integer}}) = fill_reshape(A, dims...)
reshape(A::AbstractFill, dims::Tuple{Integer,Infinity,Vararg{Integer}}) = fill_reshape(A, dims...)


BroadcastStyle(::Type{<:ReshapedArray{T,N,<:Any,NTuple{N,<:Infinity}}}) where {T,N} = LazyArrayStyle{N}()
BroadcastStyle(::Type{<:ReshapedArray{T,2,<:Any,<:Tuple{<:Any,<:Infinity}}}) where {T} = LazyArrayStyle{2}()
BroadcastStyle(::Type{<:ReshapedArray{T,2,<:Any,<:Tuple{<:Infinity,<:Any}}}) where {T} = LazyArrayStyle{2}()


MemoryLayout(::Type{<:ReshapedArray{T,N,A,DIMS}}) where {T,N,A,DIMS} = reshapedlayout(MemoryLayout(A), DIMS)


###
# permutedims for reshaped unrolls
###

permutedims(R::ReshapedArray{<:Any,2,<:AbstractVector}) = parent(R)