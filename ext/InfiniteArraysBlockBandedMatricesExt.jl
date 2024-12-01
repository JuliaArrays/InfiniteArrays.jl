module InfiniteArraysBlockBandedMatricesExt

using InfiniteArrays, BlockBandedMatrices

using InfiniteArrays.LinearAlgebra, BlockBandedMatrices.BlockArrays, BlockBandedMatrices.BandedMatrices

import BlockBandedMatrices: BlockSkylineSizes, BlockBandedMatrix, _BlockSkylineMatrix
import InfiniteArrays: InfiniteCardinal, OneToInf, InfStepRange
import BlockBandedMatrices.BandedMatrices: _BandedMatrix

const BlockTriPertToeplitz{T} = BlockMatrix{T,Tridiagonal{Matrix{T},Vcat{Matrix{T},1,Tuple{Vector{Matrix{T}},Fill{Matrix{T},1,Tuple{OneToInf{Int}}}}}},
                                        NTuple{2,BlockedOneTo{Int,Vcat{Int,1,Tuple{Vector{Int},InfStepRange{Int,Int}}}}}}


BlockBandedMatrices.blockbanded_colstop(A, x::InfiniteCardinal{0}) = x
BlockBandedMatrices.blockbanded_rowstop(A, x::InfiniteCardinal{0}) = x

function BlockSkylineSizes(A::BlockTriPertToeplitz, (l,u)::NTuple{2,Int})
    N = max(length(A.blocks.du.args[1])+1,length(A.blocks.d.args[1]),length(A.blocks.dl.args[1]))
    block_sizes = Vector{Int}(undef, N) # assume square
    block_starts = BandedMatrix{Int}(undef, (N+l,N),  (l,u))
    block_strides = Vector{Int}(undef, N)
    for J=1:N
        block_starts[max(1,J-u),J] = J == 1 ? 1 :
                            block_starts[max(1,J-1-u),J-1]+block_sizes[J-1]*block_strides[J-1]

        for K=max(1,J-u)+1:J+l
            block_starts[K,J] = block_starts[K-1,J]+size(A[Block(K-1,J)],1)
        end
        block_strides[J] = block_starts[J+l,J] + size(A[Block(J+l,J)],1) - block_starts[max(1,J-u),J]
        block_sizes[J] = size(A[Block(J,J)],2)
    end

    block_stride∞ = 0
    for K=max(1,N+1-u):N+1+l
        block_stride∞ += size(A[Block(K,N+1)],1)
    end
    block_size∞ = size(A[Block(N+1,N+1)],2)

    bs∞ = fill(block_starts[max(1,N-u),N]+block_strides[N]*size(A[Block(N,N)],2):block_stride∞*block_size∞:∞, l+u+1)
    for k=2:l+u+1
        bs∞[k] = bs∞[k-1] .+ size(A[Block(N+1-u+k-1,N+1)],1)
    end

    BlockSkylineSizes(axes(A),
                        _BandedMatrix(Hcat(block_starts.data, Vcat(adjoint.(bs∞)...)), ℵ₀, l, u),
                        Vcat(block_strides, Fill(block_stride∞,∞)),
                        Fill(l,∞),Fill(u,∞))
end

function BlockBandedMatrix(A::BlockTriPertToeplitz{T}, (l,u)::NTuple{2,Int}) where T
    data = T[]
    append!(data,vec(A[Block.(1:1+l),Block(1)]))
    N = max(length(A.blocks.du.args[1])+1,length(A.blocks.d.args[1]),length(A.blocks.dl.args[1]))
    for J=2:N
        append!(data, vec(A[Block.(max(1,J-u):J+l),Block(J)]))
    end
    tl = mortar(Fill(vec(A[Block.(max(1,N+1-u):N+1+l),Block(N+1)]),∞))

    B = _BlockSkylineMatrix(Vcat(data,tl), BlockSkylineSizes(A, (l,u)))
end

BlockBandedMatrix(A::BlockTriPertToeplitz) = BlockBandedMatrix(A, blockbandwidths(A))

end #module