struct LazyArrayStyle{N} <: AbstractArrayStyle{N} end

BroadcastStyle(::Type{<:AbstractInfUnitRange}) = LazyArrayStyle{1}()
BroadcastStyle(::Type{<:Diagonal{<:Any,<:AbstractInfUnitRange}}) = LazyArrayStyle{2}()
BroadcastStyle(::Type{<:Ones{T,N,NTuple{N,Infinity}}}) where {T,N} = LazyArrayStyle{N}()
BroadcastStyle(::Type{<:Ones{T,2,Tuple{Int,Infinity}}}) where {T} = LazyArrayStyle{2}()
BroadcastStyle(::Type{<:Ones{T,2,Tuple{Infinity,Int}}}) where {T} = LazyArrayStyle{2}()

const InfIndexRanges{T<:Integer} = Union{InfStepRange{T},AbstractInfUnitRange{T},Slice{OneToInf{T}}}

BroadcastStyle(::Type{<:SubArray{<:Any,1,<:Any,Tuple{<:InfIndexRanges}}})= LazyArrayStyle{1}()
BroadcastStyle(::Type{<:SubArray{<:Any,2,<:Any,<:Tuple{<:InfIndexRanges,<:InfIndexRanges}}})= LazyArrayStyle{1}()
BroadcastStyle(::Type{<:SubArray{<:Any,2,<:Any,<:Tuple{<:InfIndexRanges,<:Any}}})= LazyArrayStyle{1}()
BroadcastStyle(::Type{<:SubArray{<:Any,2,<:Any,<:Tuple{<:Any,<:InfIndexRanges}}})= LazyArrayStyle{1}()

struct BroadcastArray{T, N, FF, AA, INDS} <: AbstractArray{T, N}
    f::FF
    terms::AA
    axes::INDS
end

BroadcastArray(f, ::Type{ElType}, inds::Indices{N}, As...) where {ElType,N} =
    BroadcastArray{ElType,N,typeof(f),typeof(As),typeof(inds)}(f, As, inds)
BroadcastArray(f, A, Bs...) = BroadcastArray(f, combine_eltypes(f, A, Bs...), combine_indices(A, Bs...), A, Bs...)

axes(A::BroadcastArray) = A.axes
size(A::BroadcastArray) = length.(A.axes)


@generated function getindex(A::BroadcastArray{<:Any,N,<:Any,<:Tuple{Vararg{Any,M}}}, kj::CartesianIndex{N}) where {N,M}
    quote
        $(Expr(:meta, :inline))
        @nexprs $M i->(A_i = A.terms[i])
        @nexprs $M i->(@inbounds val_i = _broadcast_getindex(A_i, kj))
        # call the function and store the result
        @ncall $M A.f val
    end
end

getindex(A::BroadcastArray{<:Any,N}, kj::Vararg{Integer,N}) where N = getindex(A, CartesianIndex(kj...))


@inline broadcast(f, s::LazyArrayStyle, ::Type{ElType}, inds::Indices, As...) where ElType =
    BroadcastArray(f, ElType, inds, As...)
