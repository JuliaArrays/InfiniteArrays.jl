## CachedOperator

struct CachedArray{T,N,DM<:AbstractArray{T,N},M<:AbstractArray{T,N}} <: AbstractArray{T,N}
    data::DM
    array::M
    datasize::NTuple{N,Int}
end

const CachedVector{T,DM<:AbstractVector{T},M<:AbstractVector{T}} = CachedArray{T,1,DM,M}
const CachedMatrix{T,DM<:AbstractMatrix{T},M<:AbstractMatrix{T}} = CachedArray{T,2,DM,M}

CachedArray(data::AbstractArray{T,N}, array::AbstractArray{T,N}, sz::NTuple{N,Int}) where {T,N} =
    CachedArray{T,N,typeof(data),typeof(array)}(array, data, sz)
CachedArray(data::AbstractArray, array::AbstractArray) = CachedArray(array, data, size(data))

CachedArray(::Type{Array}, array::AbstractArray{T,N}) where {T,N} =
    CachedArray(Array{T,N}(undef, ntuple(zero,N)), array)


CachedArray(array::AbstractArray{T,N}) where {T,N} =
    CachedArray(similar(array, ntuple(zero,N)), array)

doc"""
    cache(array::AbstractArray)

Caches the entries of an array.
"""
cache(O::AbstractArray) = CachedArray(O)
cache(::Type{MT}, O::AbstractArray) where {MT<:AbstractArray} = CachedArray(MT,O;kwds...)

convert(::Type{AbstractArray{T}}, S::CachedArray{T}) where {T} = S
convert(::Type{AbstractArray{T}}, S::CachedArray) where {T} =
    CachedArray(convert(AbstractArray{T}, S.data), convert(AbstractArray{T}, S.array), S.datasize)


size(A::CachedArray) = size(A.array)

@propagate_inbounds function Base.getindex(B::CachedArray{T,N}, kj::Vararg{Integer,N}) where {T,N}
    @boundscheck checkbounds(Bool, B, kj)
    resizedata!(B, kj...)
    B.data[kj...]
end

@propagate_inbounds function Base.setindex!(B::CachedArray{T,N}, v, kj::Vararg{Integer,N}) where {T,N}
    @boundscheck checkbounds(Bool, B, kj)
    resizedata!(B,kj...)
    @inbounds B.data[kj...] = v
    v
end


## Array caching

function resizedata!(B::CachedOperator{T,N,Array{T,N}},nm::Vararg{Integer,N}) where {T<:Number,N}
    if n > size(B,1) || m > size(B,2)
        throw(ArgumentError("Cannot resize beyound size of operator"))
    end

    # this does nothing if already in dimensions
    N,M=size(B.data)
    if n > N && m > M
        B.data = unsafe_resize!(B.data,n,m)
    elseif n > N
        B.data = unsafe_resize!(B.data,n,:)
    elseif m > M
        B.data = unsafe_resize!(B.data,:,m)
    end

    if n ≤ B.datasize[1] && m ≤ B.datasize[2]
        # do nothing
        B
    elseif n ≤ B.datasize[1]
        kr,jr=1:B.datasize[1],B.datasize[2]+1:m
        B.data[kr,jr] = B.op[kr,jr]
        B.datasize = (B.datasize[1],m)
        B
    elseif m ≤ B.datasize[2]
        kr,jr=B.datasize[1]+1:n,1:B.datasize[2]
        B.data[kr,jr] = B.op[kr,jr]
        B.datasize = (n,B.datasize[2])
        B
    else
        # resize rows then columns
        resizedata!(resizedata!(B,n,B.datasize[2]),n,m)
    end
end
