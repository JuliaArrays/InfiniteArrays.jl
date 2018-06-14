struct LazyArrayStyle{N} <: AbstractArrayStyle{N} end

BroadcastStyle(::Type{<:AbstractInfUnitRange}) = LazyArrayStyle{1}()
BroadcastStyle(::Type{<:Diagonal{<:Any,<:AbstractInfUnitRange}}) = LazyArrayStyle{2}()
for typ in (:Ones, :Zeros, :Fill)
    @eval begin
        BroadcastStyle(::Type{$typ{T,N,NTuple{N,Infinity}}}) where {T,N} = LazyArrayStyle{N}()
        BroadcastStyle(::Type{$typ{T,2,Tuple{Int,Infinity}}}) where {T} = LazyArrayStyle{2}()
        BroadcastStyle(::Type{$typ{T,2,Tuple{Infinity,Int}}}) where {T} = LazyArrayStyle{2}()
    end
end

BroadcastStyle(::Type{Eye{T,NTuple{2,Infinity}}}) where {T} = LazyArrayStyle{2}()
BroadcastStyle(::Type{Eye{T,Tuple{Int,Infinity}}}) where {T} = LazyArrayStyle{2}()
BroadcastStyle(::Type{Eye{T,Tuple{Infinity,Int}}}) where {T} = LazyArrayStyle{2}()

const InfIndexRanges{T<:Integer} = Union{InfStepRange{T},AbstractInfUnitRange{T},Slice{OneToInf{T}}}

BroadcastStyle(::Type{<:SubArray{<:Any,1,<:Any,Tuple{<:InfIndexRanges}}})= LazyArrayStyle{1}()
BroadcastStyle(::Type{<:SubArray{<:Any,2,<:Any,<:Tuple{<:InfIndexRanges,<:InfIndexRanges}}})= LazyArrayStyle{1}()
BroadcastStyle(::Type{<:SubArray{<:Any,2,<:Any,<:Tuple{<:InfIndexRanges,<:Any}}})= LazyArrayStyle{1}()
BroadcastStyle(::Type{<:SubArray{<:Any,2,<:Any,<:Tuple{<:Any,<:InfIndexRanges}}})= LazyArrayStyle{1}()

struct BroadcastArray{T, N, BRD<:Broadcasted} <: AbstractArray{T, N}
    broadcasted::BRD
end

axes(A::BroadcastArray) = A.broadcasted.axes
size(A::BroadcastArray) = length.(A.axes)


@propagate_inbounds getindex(A::BroadcastArray{<:Any,N,<:Any,<:Tuple{Vararg{Any,M}}}, kj::CartesianIndex{N}) where {N,M} =
    A.broadcasted[kj]

getindex(A::BroadcastArray{<:Any,N}, kj::Vararg{Integer,N}) where N = getindex(A, CartesianIndex(kj...))


copy(bc::Broadcasted{<:LazyArrayStyle}) = BroadcastArray(bc)
