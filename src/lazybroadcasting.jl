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

BroadcastArray{T,N}(bc::BRD) where {T,N,BRD<:Broadcasted} = BroadcastArray{T,N,BRD}(bc)
BroadcastArray{T}(bc::Broadcasted{<:Any,<:Tuple{Vararg{Any,N}}}) where {T,N} = BroadcastArray{T,N}(bc)
BroadcastArray(bc::Broadcasted) = BroadcastArray{combine_eltypes(bc.f, bc.args)}(bc)
BroadcastArray(b::Broadcasted{<:Any,Nothing}) = BroadcastArray(instantiate(b))

axes(A::BroadcastArray) = axes(A.broadcasted)
size(A::BroadcastArray) = map(length, axes(A))


@propagate_inbounds getindex(A::BroadcastArray, kj...) = A.broadcasted[kj...]

copy(bc::Broadcasted{<:LazyArrayStyle}) = BroadcastArray(bc)
