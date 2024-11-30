module InfiniteArraysBlockArraysExt
using InfiniteArrays, BlockArrays
using InfiniteArrays.ArrayLayouts, InfiniteArrays.LazyArrays, InfiniteArrays.LinearAlgebra

import Base: length, size, axes, BroadcastStyle
import Base.Broadcast: Broadcasted
import ArrayLayouts: sub_materialize, axes_print_matrix_row
import InfiniteArrays: OneToInf, PosInfinity, InfRanges, RealInfinity, Infinity, InfStepRange
import BlockArrays: AbstractBlockLayout, sizes_from_blocks, BlockTridiagonal, OneToCumsum, BlockSlice, AbstractBlockedUnitRange,
                    BlockLayout
import LazyArrays: PaddedColumns

const OneToInfCumsum = RangeCumsum{Int,OneToInf{Int}}

BlockArrays.sortedunion(::AbstractVector{<:PosInfinity}, ::AbstractVector{<:PosInfinity}) = [∞]
function BlockArrays.sortedunion(::AbstractVector{<:PosInfinity}, b)
    @assert isinf(length(b))
    b
end

function BlockArrays.sortedunion(b, ::AbstractVector{<:PosInfinity})
    @assert isinf(length(b))
    b
end
BlockArrays.sortedunion(a::OneToInfCumsum, ::OneToInfCumsum) = a

BlockArrays.blocklasts(a::InfRanges) = Fill(length(a),1)

BlockArrays.findblock(::BlockedOneTo, ::RealInfinity) = Block(ℵ₀)

function BlockArrays.sortedunion(a::Vcat{Int,1,<:Tuple{Union{Int,AbstractVector{Int}},<:AbstractRange}},
                                 b::Vcat{Int,1,<:Tuple{Union{Int,AbstractVector{Int}},<:AbstractRange}})
    @assert a == b # TODO: generailse? Not sure how to do so in a type stable fashion
    a
end

sizes_from_blocks(A::AbstractVector, ::Tuple{OneToInf{Int}}) = (map(length,A),)
length(::BlockedOneTo{Int,<:InfRanges}) = ℵ₀

const OneToInfBlocks = BlockedOneTo{Int,OneToInfCumsum}
const OneToBlocks = BlockedOneTo{Int,OneToCumsum}

axes(a::OneToInfBlocks) = (a,)
axes(a::OneToBlocks) = (a,)


sub_materialize(_, V, ::Tuple{BlockedOneTo{Int,<:InfRanges}}) = V
sub_materialize(::AbstractBlockLayout, V, ::Tuple{BlockedOneTo{Int,<:InfRanges}}) = V
function sub_materialize(::PaddedColumns, v::AbstractVector{T}, ax::Tuple{BlockedOneTo{Int,<:InfRanges}}) where T
    dat = paddeddata(v)
    BlockedVector(Vcat(sub_materialize(dat), Zeros{T}(∞)), ax)
end

BlockArrays.dimlength(start, ::Infinity) = ℵ₀

function copy(bc::Broadcasted{<:BroadcastStyle,<:Any,typeof(*),<:Tuple{Ones{T,1,Tuple{OneToInfBlocks}},AbstractArray{V,N}}}) where {N,T,V}
    a,b = bc.args
    @assert bc.axes == axes(b)
    convert(AbstractArray{promote_type(T,V),N}, b)
end

function copy(bc::Broadcasted{<:BroadcastStyle,<:Any,typeof(*),<:Tuple{AbstractArray{T,N},Ones{V,1,Tuple{OneToInfBlocks}}}}) where {N,T,V}
    a,b = bc.args
    @assert bc.axes == axes(a)
    convert(AbstractArray{promote_type(T,V),N}, a)
end

_block_interlace_axes(::Int, ax::Tuple{BlockedOneTo{Int,OneToInf{Int}}}...) = (blockedrange(Fill(length(ax), ∞)),)

_block_interlace_axes(nbc::Int, ax::NTuple{2,BlockedOneTo{Int,OneToInf{Int}}}...) =
    (blockedrange(Fill(length(ax) ÷ nbc, ∞)),blockedrange(Fill(mod1(length(ax),nbc), ∞)))

#######
# block broadcasted
######


BroadcastStyle(::Type{<:SubArray{T,N,Arr,<:NTuple{N,BlockSlice{BlockRange{1,Tuple{II}}}},false}}) where {T,N,Arr<:BlockArray,II<:InfRanges} =
    LazyArrayStyle{N}()

# TODO: generalise following
BroadcastStyle(::Type{<:BlockArray{T,N,<:AbstractArray{<:AbstractArray{T,N},N},<:NTuple{N,BlockedOneTo{Int,<:InfRanges}}}}) where {T,N} = LazyArrayStyle{N}()
# BroadcastStyle(::Type{<:BlockedArray{T,N,<:AbstractArray{T,N},<:NTuple{N,BlockedOneTo{Int,<:InfRanges}}}}) where {T,N} = LazyArrayStyle{N}()
BroadcastStyle(::Type{<:BlockArray{T,N,<:AbstractArray{<:AbstractArray{T,N},N},<:NTuple{N,BlockedOneTo{Int,<:RangeCumsum{Int,<:InfRanges}}}}}) where {T,N} = LazyArrayStyle{N}()
# BroadcastStyle(::Type{<:BlockedArray{T,N,<:AbstractArray{T,N},<:NTuple{N,BlockedOneTo{Int,<:RangeCumsum{Int,<:InfRanges}}}}}) where {T,N} = LazyArrayStyle{N}()


# Block banded support

sizes_from_blocks(A::Diagonal, ::NTuple{2,OneToInf{Int}}) = size.(A.diag, 1), size.(A.diag,2)
sizes_from_blocks(A::Tridiagonal, ::NTuple{2,OneToInf{Int}}) = size.(A.d, 1), size.(A.d,2)
sizes_from_blocks(A::Bidiagonal, ::NTuple{2,OneToInf{Int}}) = size.(A.dv, 1), size.(A.dv,2)


axes_print_matrix_row(::NTuple{2,AbstractBlockedUnitRange}, io, X, A, i, ::AbstractVector{<:PosInfinity}, sep) = nothing


const BlockTriPertToeplitz{T} = BlockMatrix{T,Tridiagonal{Matrix{T},Vcat{Matrix{T},1,Tuple{Vector{Matrix{T}},Fill{Matrix{T},1,Tuple{OneToInf{Int}}}}}},
                                        NTuple{2,BlockedOneTo{Int,Vcat{Int,1,Tuple{Vector{Int},InfStepRange{Int,Int}}}}}}

const BlockTridiagonalToeplitzLayout = BlockLayout{TridiagonalToeplitzLayout,DenseColumnMajor}

function BlockTridiagonal(adjA::Adjoint{T,BlockTriPertToeplitz{T}}) where T
    A = parent(adjA)
    BlockTridiagonal(Matrix.(adjoint.(A.blocks.du)),
                     Matrix.(adjoint.(A.blocks.d)),
                     Matrix.(adjoint.(A.blocks.dl)))
end

for op in (:-, :+)
    @eval begin
        function $op(A::BlockTriPertToeplitz{T}, λ::UniformScaling) where T
            TV = promote_type(T,eltype(λ))
            BlockTridiagonal(Vcat(convert.(AbstractVector{Matrix{TV}}, A.blocks.dl.args)...),
                             Vcat(convert.(AbstractVector{Matrix{TV}}, broadcast($op, A.blocks.d, Ref(λ)).args)...),
                             Vcat(convert.(AbstractVector{Matrix{TV}}, A.blocks.du.args)...))
        end
        function $op(λ::UniformScaling, A::BlockTriPertToeplitz{V}) where V
            TV = promote_type(eltype(λ),V)
            BlockTridiagonal(Vcat(convert.(AbstractVector{Matrix{TV}}, broadcast($op, A.blocks.dl.args))...),
                             Vcat(convert.(AbstractVector{Matrix{TV}}, broadcast($op, Ref(λ), A.blocks.d).args)...),
                             Vcat(convert.(AbstractVector{Matrix{TV}}, broadcast($op, A.blocks.du.args))...))
        end
        $op(adjA::Adjoint{T,BlockTriPertToeplitz{T}}, λ::UniformScaling) where T = $op(BlockTridiagonal(adjA), λ)
        $op(λ::UniformScaling, adjA::Adjoint{T,BlockTriPertToeplitz{T}}) where T = $op(λ, BlockTridiagonal(adjA))
    end
end


end # module