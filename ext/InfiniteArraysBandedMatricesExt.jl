module InfiniteArraysBandedMatricesExt
using InfiniteArrays, BandedMatrices, LinearAlgebra
using InfiniteArrays.LazyArrays, InfiniteArrays.ArrayLayouts, InfiniteArrays.FillArrays

import Base: BroadcastStyle, size, getindex, similar, copy, *, +, -, /, \, materialize!, copyto!, OneTo
import Base.Broadcast: Broadcasted
import InfiniteArrays: InfIndexRanges, Infinity, PosInfinity, OneToInf, InfAxes, AbstractInfUnitRange, InfRanges
import ArrayLayouts: sub_materialize, MemoryLayout, sublayout, mulreduce, triangularlayout, MatLdivVec, subdiagonaldata, diagonaldata, supdiagonaldata
import LazyArrays: applybroadcaststyle, applylayout, islazy, islazy_layout, simplifiable, AbstractLazyLayout, PaddedColumns, LazyArrayStyle, ApplyLayout, AbstractLazyBandedLayout, ApplyBandedLayout, BroadcastBandedLayout
import BandedMatrices: _BandedMatrix, AbstractBandedMatrix, banded_similar, BandedMatrix, bandedcolumns, BandedColumns, bandeddata
import FillArrays: AbstractFillMatrix, AbstractFill, getindex_value

BroadcastStyle(::Type{<:SubArray{<:Any,2,<:AbstractBandedMatrix,<:Tuple{<:InfIndexRanges,<:InfIndexRanges}}})= LazyArrayStyle{2}()
_BandedMatrix(data::AbstractMatrix{T}, ::Infinity, l, u) where T = _BandedMatrix(data, ℵ₀, l, u)

###
# BandIndexing
###

struct InfBandCartesianIndices <: AbstractVector{CartesianIndex{2}}
    b::Int
end

InfBandCartesianIndices(b::Band) = InfBandCartesianIndices(b.i)

size(::InfBandCartesianIndices) = (∞,)
getindex(B::InfBandCartesianIndices, k::Int) = B.b ≥ 0 ? CartesianIndex(k, k+B.b) : CartesianIndex(k-B.b, k)

Base.checkindex(::Type{Bool}, ::NTuple{2,OneToInf{Int}}, ::InfBandCartesianIndices) = true
BandedMatrices.band_to_indices(_, ::NTuple{2,OneToInf{Int}}, b) = (InfBandCartesianIndices(b),)
BroadcastStyle(::Type{<:SubArray{<:Any,1,<:Any,Tuple{InfBandCartesianIndices}}}) = LazyArrayStyle{1}()

_inf_banded_sub_materialize(_, V) = V
function _inf_banded_sub_materialize(::BandedColumns, V)
    A = parent(V)
    b = parentindices(V)[1].b
    data = bandeddata(A)
    l,u = bandwidths(A)
    if -l ≤ b ≤ u
        data[u+1-b, max(1,b+1):end]
    else
        Zeros{eltype(V)}(∞) # Not type stable
    end
end

sub_materialize(_, V::SubArray{<:Any,1,<:AbstractMatrix,Tuple{InfBandCartesianIndices}}, ::Tuple{InfAxes}) =
    _inf_banded_sub_materialize(MemoryLayout(parent(V)), V)







const TriToeplitz{T} = Tridiagonal{T,Fill{T,1,Tuple{OneToInf{Int}}}}
const ConstRowMatrix{T} = ApplyMatrix{T,typeof(*),<:Tuple{<:AbstractVector,<:AbstractFillMatrix{<:Any,Tuple{OneTo{Int},OneToInf{Int}}}}}
const PertConstRowMatrix{T} = Hcat{T,<:Tuple{Array{T},<:ConstRowMatrix{T}}}
const InfToeplitz{T} = BandedMatrix{T,<:ConstRowMatrix{T},OneToInf{Int}}
const PertToeplitz{T} = BandedMatrix{T,<:PertConstRowMatrix{T},OneToInf{Int}}

const SymTriPertToeplitz{T} = SymTridiagonal{T,Vcat{T,1,Tuple{Vector{T},Fill{T,1,Tuple{OneToInf{Int}}}}}}
const TriPertToeplitz{T} = Tridiagonal{T,Vcat{T,1,Tuple{Vector{T},Fill{T,1,Tuple{OneToInf{Int}}}}}}
const AdjTriPertToeplitz{T} = Adjoint{T,Tridiagonal{T,Vcat{T,1,Tuple{Vector{T},Fill{T,1,Tuple{OneToInf{Int}}}}}}}
const InfBandedMatrix{T,V<:AbstractMatrix{T}} = BandedMatrix{T,V,OneToInf{Int}}

_prepad(p, a) = Vcat(Zeros{eltype(a)}(max(p,0)), a)
_prepad(p, a::Zeros{T,1}) where T = Zeros{T}(length(a)+p)
_prepad(p, a::Ones{T,1}) where T = Ones{T}(length(a)+p)
_prepad(p, a::AbstractFill{T,1}) where T = Fill{T}(getindex_value(a), length(a)+p)

banded_similar(T, (m,n)::Tuple{Int,PosInfinity}, (l,u)::Tuple{Int,Int}) = BandedMatrix{T}(undef, (n,m), (u,l))'

function BandedMatrix{T}(kv::Tuple{Vararg{Pair{<:Integer,<:AbstractVector}}},
                         ::NTuple{2,PosInfinity},
                         (l,u)::NTuple{2,Integer}) where T
    ks = getproperty.(kv, :first)
    rws = Vcat(permutedims.(_prepad.(ks,getproperty.(kv, :second)))...)
    c = zeros(T, l+u+1, length(ks))
    for (k,j) in zip(u .- ks .+ 1,1:length(ks))
        c[k,j] = one(T)
    end
    _BandedMatrix(ApplyArray(*,c,rws), ℵ₀, l, u)
end

# Construct InfToeplitz
function BandedMatrix{T}(kv::Tuple{Vararg{Pair{<:Integer,<:Fill{<:Any,1,Tuple{OneToInf{Int}}}}}},
                         mn::NTuple{2,PosInfinity},
                         lu::NTuple{2,Integer}) where T
    m,n = mn
    @assert isinf(n)
    l,u = lu
    t = zeros(T, u+l+1)
    for (k,v) in kv
        p = length(v)
        t[u-k+1] = v.value
    end

    return _BandedMatrix(t * Ones{T}(1,∞), Integer(m), l, u)
end

function BandedMatrix{T}(kv::Tuple{Vararg{Pair{<:Integer,<:Vcat{<:Any,1,<:Tuple{<:AbstractVector,Fill{<:Any,1,Tuple{OneToInf{Int}}}}}}}},
                         mn::NTuple{2,PosInfinity},
                         lu::NTuple{2,Integer}) where T
    m,n = mn
    @assert isinf(n)
    l,u = lu
    M = mapreduce(x -> length(x.second.args[1]) + max(0,x.first), max, kv) # number of data rows
    data = zeros(T, u+l+1, M)
    t = zeros(T, u+l+1)
    for (k,v) in kv
        a,b = v.args
        p = length(a)
        t[u-k+1] = b.value
        if k ≤ 0
            data[u-k+1,1:p] = a
            data[u-k+1,p+1:end] .= b.value
        else
            data[u-k+1,k+1:k+p] = a
            data[u-k+1,k+p+1:end] .= b.value
        end
    end

    return _BandedMatrix(Hcat(data, t * Ones{T}(1,∞)), Integer(m), l, u)
end


function BandedMatrix(Ac::Adjoint{T,<:InfToeplitz}) where T
    A = parent(Ac)
    l,u = bandwidths(A)
    a = A.data.args[1]
    _BandedMatrix(reverse(conj(a)) * Ones{T}(1,∞), ℵ₀, u, l)
end

function BandedMatrix(Ac::Transpose{T,<:InfToeplitz}) where T
    A = parent(Ac)
    l,u = bandwidths(A)
    a = A.data.args[1]
    _BandedMatrix(reverse(a) * Ones{T}(1,∞), ℵ₀, u, l)
end

function BandedMatrix(Ac::Adjoint{T,<:PertToeplitz}) where T
    A = parent(Ac)
    l,u = bandwidths(A)
    a,b = A.data.args
    Ac_fd = BandedMatrix(_BandedMatrix(Hcat(a, b[:,1:l+1]), size(a,2)+l, l, u)')
    _BandedMatrix(Hcat(Ac_fd.data, reverse(conj(b.args[1])) * Ones{T}(1,∞)), ℵ₀, u, l)
end

function BandedMatrix(Ac::Transpose{T,<:PertToeplitz}) where T
    A = parent(Ac)
    l,u = bandwidths(A)
    a,b = A.data.args
    Ac_fd = BandedMatrix(transpose(_BandedMatrix(Hcat(a, b[:,1:l+1]), size(a,2)+l, l, u)))
    _BandedMatrix(Hcat(Ac_fd.data, reverse(b.args[1]) * Ones{T}(1,∞)), ℵ₀, u, l)
end


for op in (:-, :+)
    @eval begin
        function $op(A::SymTriPertToeplitz{T}, λ::UniformScaling) where T
            TV = promote_type(T,eltype(λ))
            dv = Vcat(convert.(AbstractVector{TV}, A.dv.args)...)
            ev = Vcat(convert.(AbstractVector{TV}, A.ev.args)...)
            SymTridiagonal(broadcast($op, dv, Ref(λ.λ)), ev)
        end
        function $op(λ::UniformScaling, A::SymTriPertToeplitz{V}) where V
            TV = promote_type(eltype(λ),V)
            SymTridiagonal(convert(AbstractVector{TV}, broadcast($op, Ref(λ.λ), A.dv)),
                           convert(AbstractVector{TV}, broadcast($op, A.ev)))
        end
        function $op(A::SymTridiagonal{T,<:AbstractFill}, λ::UniformScaling) where T
            TV = promote_type(T,eltype(λ))
            SymTridiagonal(convert(AbstractVector{TV}, broadcast($op, A.dv, Ref(λ.λ))),
                           convert(AbstractVector{TV}, A.ev))
        end

        function $op(A::TriPertToeplitz{T}, λ::UniformScaling) where T
            TV = promote_type(T,eltype(λ))
            Tridiagonal(Vcat(convert.(AbstractVector{TV}, A.dl.args)...),
                        Vcat(convert.(AbstractVector{TV}, broadcast($op, A.d, λ.λ).args)...),
                        Vcat(convert.(AbstractVector{TV}, A.du.args)...))
        end
        function $op(λ::UniformScaling, A::TriPertToeplitz{V}) where V
            TV = promote_type(eltype(λ),V)
            Tridiagonal(Vcat(convert.(AbstractVector{TV}, broadcast($op, A.dl.args))...),
                        Vcat(convert.(AbstractVector{TV}, broadcast($op, λ.λ, A.d).args)...),
                        Vcat(convert.(AbstractVector{TV}, broadcast($op, A.du.args))...))
        end
        function $op(adjA::AdjTriPertToeplitz{T}, λ::UniformScaling) where T
            A = parent(adjA)
            TV = promote_type(T,eltype(λ))
            Tridiagonal(Vcat(convert.(AbstractVector{TV}, A.du.args)...),
                        Vcat(convert.(AbstractVector{TV}, broadcast($op, A.d, λ.λ).args)...),
                        Vcat(convert.(AbstractVector{TV}, A.dl.args)...))
        end
        function $op(λ::UniformScaling, adjA::AdjTriPertToeplitz{V}) where V
            A = parent(adjA)
            TV = promote_type(eltype(λ),V)
            Tridiagonal(Vcat(convert.(AbstractVector{TV}, broadcast($op, A.du.args))...),
                        Vcat(convert.(AbstractVector{TV}, broadcast($op, λ.λ, A.d).args)...),
                        Vcat(convert.(AbstractVector{TV}, broadcast($op, A.dl.args))...))
        end

        function $op(λ::UniformScaling, A::InfToeplitz{V}) where V
            l,u = bandwidths(A)
            TV = promote_type(eltype(λ),V)
            a = convert(AbstractVector{TV}, $op.(A.data.args[1]))
            a[u+1] += λ.λ
            _BandedMatrix(a*Ones{TV}(1,∞), ℵ₀, l, u)
        end

        function $op(A::InfToeplitz{T}, λ::UniformScaling) where T
            l,u = bandwidths(A)
            TV = promote_type(T,eltype(λ))
            a = TV[Zeros{TV}(max(-u,0)); A.data.args[1]; Zeros{TV}(max(-l,0))]
            a[max(0,u)+1] = $op(a[max(u,0)+1], λ.λ)
            _BandedMatrix(a*Ones{TV}(1,∞), ℵ₀, max(l,0), max(u,0))
        end

        function $op(λ::UniformScaling, A::PertToeplitz{V}) where V
            l,u = bandwidths(A)
            TV = promote_type(eltype(λ),V)
            a, t = map(AbstractArray{TV}, A.data.args)
            b = $op.(t.args[1])
            a[u+1,:] .+= λ.λ
            b[u+1] += λ.λ
            _BandedMatrix(Hcat(a, b*Ones{TV}(1,∞)), ℵ₀, l, u)
        end

        function $op(A::PertToeplitz{T}, λ::UniformScaling) where T
            l,u = bandwidths(A)
            TV = promote_type(T,eltype(λ))
            ã, t = A.data.args
            a = AbstractArray{TV}(ã)
            b = AbstractVector{TV}(t.args[1])
            a[u+1,:] .= $op.(a[u+1,:],λ.λ)
            b[u+1] = $op(b[u+1], λ.λ)
            _BandedMatrix(Hcat(a, b*Ones{TV}(1,∞)), ℵ₀, l, u)
        end
    end
end



####
# Conversions to BandedMatrix
####

function BandedMatrix(A::PertToeplitz{T}, (l,u)::Tuple{Int,Int}) where T
    @assert A.u == u # Not implemented
    a, b = A.data.args
    t = b.args[1] # topelitz part
    t_pad = vcat(t,Zeros(l-A.l))
    data = Hcat([vcat(a,Zeros{T}(l-A.l,size(a,2))) repeat(t_pad,1,l)], t_pad * Ones{T}(1,∞))
    _BandedMatrix(data, ℵ₀, l, u)
end

function BandedMatrix(A::SymTriPertToeplitz{T}, (l,u)::Tuple{Int,Int}) where T
    a,a∞ = A.dv.args
    b,b∞ = A.ev.args
    n = max(length(a), length(b)+1) + 1
    data = zeros(T, l+u+1, n)
    data[u,2:length(b)+1] .= b
    data[u,length(b)+2:end] .= b∞.value
    data[u+1,1:length(a)] .= a
    data[u+1,length(a)+1:end] .= a∞.value
    data[u+2,1:length(b)] .= b
    data[u+2,length(b)+1:end] .= b∞.value
    _BandedMatrix(Hcat(data, [Zeros{T}(u-1); b∞.value; a∞.value; b∞.value; Zeros{T}(l-1)] * Ones{T}(1,∞)), ℵ₀, l, u)
end

function BandedMatrix(A::SymTridiagonal{T,Fill{T,1,Tuple{OneToInf{Int}}}}, (l,u)::Tuple{Int,Int}) where T
    a∞ = A.dv
    b∞ = A.ev
    n = 2
    data = zeros(T, l+u+1, n)
    data[u,2:end] .= b∞.value
    data[u+1,1:end] .= a∞.value
    data[u+2,1:end] .= b∞.value
    _BandedMatrix(Hcat(data, [Zeros{T}(u-1); b∞.value; a∞.value; b∞.value; Zeros{T}(l-1)] * Ones{T}(1,∞)), ℵ₀, l, u)
end

function BandedMatrix(A::TriPertToeplitz{T}, (l,u)::Tuple{Int,Int}) where T
    a,a∞ = A.d.args
    b,b∞ = A.du.args
    c,c∞ = A.dl.args
    n = max(length(a), length(b)+1, length(c)-1) + 1
    data = zeros(T, l+u+1, n)
    data[u,2:length(b)+1] .= b
    data[u,length(b)+2:end] .= b∞.value
    data[u+1,1:length(a)] .= a
    data[u+1,length(a)+1:end] .= a∞.value
    data[u+2,1:length(c)] .= c
    data[u+2,length(c)+1:end] .= c∞.value
    _BandedMatrix(Hcat(data, [Zeros{T}(u-1); b∞.value; a∞.value; c∞.value; Zeros{T}(l-1)] * Ones{T}(1,∞)), ℵ₀, l, u)
end

function BandedMatrix(A::Tridiagonal{T,Fill{T,1,Tuple{OneToInf{Int}}}}, (l,u)::Tuple{Int,Int}) where T
    a∞ = A.d
    b∞ = A.du
    c∞ = A.dl
    n = 2
    data = zeros(T, l+u+1, n)
    data[u,2:end] .= b∞.value
    data[u+1,1:end] .= a∞.value
    data[u+2,1:end] .= c∞.value
    _BandedMatrix(Hcat(data, [Zeros{T}(u-1); b∞.value; a∞.value; c∞.value; Zeros{T}(l-1)] * Ones{T}(1,∞)), ℵ₀, l, u)
end

function InfToeplitz(A::Tridiagonal{T,Fill{T,1,Tuple{OneToInf{Int}}}}, (l,u)::Tuple{Int,Int}) where T
    a∞ = A.d
    b∞ = A.du
    c∞ = A.dl
    _BandedMatrix([Zeros{T}(u-1); b∞.value; a∞.value; c∞.value; Zeros{T}(l-1)] * Ones{T}(1,∞), ℵ₀, l, u)
end

InfToeplitz(A::Tridiagonal{T,Fill{T,1,Tuple{OneToInf{Int}}}}) where T = InfToeplitz(A, bandwidths(A))


####
# Toeplitz layout
####

_pertdata(A::ConstRowMatrix{T}) where T = Array{T}(undef,size(A,1),0)
_pertdata(A::Hcat{T,<:Tuple{Vector{T},<:ConstRowMatrix{T}}}) where T = 1
_pertdata(A::PertConstRowMatrix) = A.args[1]
function _pertdata(A::SubArray)
    P = parent(A)
    kr,jr = parentindices(A)
    dat = _pertdata(P)
    dat[kr,jr ∩ axes(dat,2)]
end

_constrows(A::ConstRowMatrix) = A.args[1]*getindex_value(A.args[2])
_constrows(A::PertConstRowMatrix) = _constrows(A.args[2])
_constrows(A::SubArray) = _constrows(parent(A))[parentindices(A)[1]]

ConstRowMatrix(A::AbstractMatrix{T}) where T = ApplyMatrix(*, A[:,1], Ones{T}(1,size(A,2)))
PertConstRowMatrix(A::AbstractMatrix{T}) where T =
    Hcat(_pertdata(A), ApplyMatrix(*, _constrows(A), Ones{T}(1,size(A,2))))

struct ConstRows <: AbstractLazyLayout end
struct PertConstRows <: AbstractLazyLayout end
MemoryLayout(::Type{<:ConstRowMatrix}) = ConstRows()
MemoryLayout(::Type{<:PertConstRowMatrix}) = PertConstRows()
bandedcolumns(::ConstRows) = BandedToeplitzLayout()
bandedcolumns(::PertConstRows) = PertToeplitzLayout()
sublayout(::ConstRows, inds...) = sublayout(ApplyLayout{typeof(*)}(), inds...)
sublayout(::PertConstRows, inds...) = sublayout(ApplyLayout{typeof(hcat)}(), inds...)
for Typ in (:ConstRows, :PertConstRows)
    @eval begin
        sublayout(::$Typ, ::Type{<:Tuple{Any,AbstractInfUnitRange{Int}}}) = $Typ() # no way to lose const rows
        applybroadcaststyle(::Type{<:AbstractMatrix}, ::$Typ) = LazyArrayStyle{2}()
        applylayout(::Type, ::$Typ, _...) = LazyLayout()
    end
end

"""
    TridiagonalToeplitzLayout

represents a matrix which is tridiagonal and toeplitz. Must support
`subdiagonalconstant`, `diagonalconstant`, `supdiagonalconstant`.
"""
struct TridiagonalToeplitzLayout <: AbstractLazyBandedLayout end
const BandedToeplitzLayout = BandedColumns{ConstRows}
const PertToeplitzLayout = BandedColumns{PertConstRows}
const PertTriangularToeplitzLayout{UPLO,UNIT} = TriangularLayout{UPLO,UNIT,BandedColumns{PertConstRows}}
struct BidiagonalToeplitzLayout <: AbstractLazyBandedLayout end
struct PertBidiagonalToeplitzLayout <: AbstractLazyBandedLayout end
struct PertTridiagonalToeplitzLayout <: AbstractLazyBandedLayout end

const InfToeplitzLayouts = Union{TridiagonalToeplitzLayout, BandedToeplitzLayout, BidiagonalToeplitzLayout,
                                 PertToeplitzLayout, PertTriangularToeplitzLayout, PertBidiagonalToeplitzLayout, PertTridiagonalToeplitzLayout}

subdiagonalconstant(A) = getindex_value(subdiagonaldata(A))
diagonalconstant(A) = getindex_value(diagonaldata(A))
supdiagonalconstant(A) = getindex_value(supdiagonaldata(A))


islazy_layout(::InfToeplitzLayouts) = Val(true)
islazy(::BandedMatrix{<:Any,<:Any,OneToInf{Int}}) = Val(true)



_BandedMatrix(::BandedToeplitzLayout, A::AbstractMatrix) =
    _BandedMatrix(ConstRowMatrix(bandeddata(A)), size(A,1), bandwidths(A)...)
_BandedMatrix(::PertToeplitzLayout, A::AbstractMatrix) =
    _BandedMatrix(PertConstRowMatrix(bandeddata(A)), size(A,1), bandwidths(A)...)

# for Lay in (:BandedToeplitzLayout, :PertToeplitzLayout)
#     @eval begin
#         sublayout(::$Lay, ::Type{<:Tuple{AbstractInfUnitRange{Int},AbstractInfUnitRange{Int}}}) = $Lay()
#         sublayout(::$Lay, ::Type{<:Tuple{Slice,AbstractInfUnitRange{Int}}}) = $Lay()
#         sublayout(::$Lay, ::Type{<:Tuple{AbstractInfUnitRange{Int},Slice}}) = $Lay()
#         sublayout(::$Lay, ::Type{<:Tuple{Slice,Slice}}) = $Lay()

#         sub_materialize(::$Lay, V) = BandedMatrix(V)
#     end
# end


@inline sub_materialize(::ApplyBandedLayout{typeof(*)}, V, ::Tuple{InfAxes,InfAxes}) = V
@inline sub_materialize(::BroadcastBandedLayout, V, ::Tuple{InfAxes,InfAxes}) = V
@inline sub_materialize(::BandedColumns, V, ::Tuple{InfAxes,InfAxes}) = BandedMatrix(V)
@inline sub_materialize(::BandedColumns, V, ::Tuple{InfAxes,OneTo{Int}}) = BandedMatrix(V)

##
# UniformScaling
##

# for op in (:+, :-), Typ in (:(BandedMatrix{<:Any,<:Any,OneToInf{Int}}),
#                             :(Adjoint{<:Any,<:BandedMatrix{<:Any,<:Any,OneToInf{Int}}}),
#                             :(Transpose{<:Any,<:BandedMatrix{<:Any,<:Any,OneToInf{Int}}}))
#     @eval begin
#         $op(A::$Typ, λ::UniformScaling) = $op(A, Diagonal(Fill(λ.λ,∞)))
#         $op(λ::UniformScaling, A::$Typ) = $op(Diagonal(Fill(λ.λ,∞)), A)
#     end
# end


_default_banded_broadcast(bc::Broadcasted, ::Tuple{<:OneToInf,<:Any}) = copy(Broadcasted{LazyArrayStyle{2}}(bc.f, bc.args))

###
# Banded * Banded
###

BandedMatrix{T}(::UndefInitializer, axes::Tuple{OneToInf{Int},OneTo{Int}}, lu::NTuple{2,Integer}) where T = BandedMatrix{T}(undef, map(length,axes), lu)
BandedMatrix{T}(::UndefInitializer, (m,n)::Tuple{Infinity,Int}, lu::NTuple{2,Integer}) where T = BandedMatrix{T}(undef, (ℵ₀,n), lu)

similar(M::MulAdd{<:AbstractBandedLayout,<:AbstractBandedLayout}, ::Type{T}, axes::Tuple{OneTo{Int},OneToInf{Int}}) where T =
    transpose(BandedMatrix{T}(undef, reverse(axes), reverse(bandwidths(M))))
similar(M::MulAdd{<:AbstractBandedLayout,<:AbstractBandedLayout}, ::Type{T}, axes::Tuple{OneToInf{Int},OneTo{Int}}) where T =
    BandedMatrix{T}(undef, axes, bandwidths(M))


###
# BandedFill * BandedFill
###

simplifiable(::Mul{<:BandedColumns{<:AbstractFillLayout},<:BandedColumns{<:AbstractFillLayout}}) = Val(true)
copy(M::MulAdd{<:BandedColumns{<:AbstractFillLayout},<:BandedColumns{<:AbstractFillLayout},ZerosLayout}) =
    _bandedfill_mul(M, axes(M.A), axes(M.B))

_bandedfill_mul(M, _, _) = copyto!(similar(M), M)
function _bandedfill_mul(M::MulAdd, ::Tuple{InfAxes,InfAxes}, ::Tuple{InfAxes,InfAxes})
    A, B = M.A, M.B
    Al,Au = bandwidths(A)
    Bl,Bu = bandwidths(B)
    l,u = Al+Bl,Au+Bu
    m = min(Au+Al,Bl+Bu)+1
    λ = getindex_value(bandeddata(A))*getindex_value(bandeddata(B))
    ret = _BandedMatrix(Hcat(Array{typeof(λ)}(undef, l+u+1,max(0,u)), [1:m-1; Fill(m,l+u-2m+3); m-1:-1:1]*Fill(λ,1,∞)), ℵ₀, l, u)
    mul!(view(ret, 1:l+u,1:u), view(A,1:l+u,1:u+Bl), view(B,1:u+Bl,1:u))
    ret
end

_bandedfill_mul(M::MulAdd, ::Tuple{InfAxes,Any}, ::Tuple{Any,Any}) = ApplyArray(*, M.A, M.B)
_bandedfill_mul(M::MulAdd, ::Tuple{InfAxes,Any}, ::Tuple{Any,InfAxes}) = ApplyArray(*, M.A, M.B)
_bandedfill_mul(M::MulAdd, ::Tuple{Any,Any}, ::Tuple{Any,InfAxes}) = ApplyArray(*, M.A, M.B)
_bandedfill_mul(M::MulAdd, ::Tuple{Any,InfAxes}, ::Tuple{InfAxes,InfAxes}) = ApplyArray(*, M.A, M.B)
_bandedfill_mul(M::MulAdd, ::Tuple{InfAxes,InfAxes}, ::Tuple{InfAxes,Any}) = ApplyArray(*, M.A, M.B)
_bandedfill_mul(M::MulAdd, ::Tuple{Any,InfAxes}, ::Tuple{InfAxes,Any}) = ApplyArray(*, M.A, M.B)

mulreduce(M::Mul{<:InfToeplitzLayouts, <:InfToeplitzLayouts}) = ApplyArray(M)
mulreduce(M::Mul{<:InfToeplitzLayouts}) = ApplyArray(M)
mulreduce(M::Mul{<:InfToeplitzLayouts,<:PaddedColumns}) = MulAdd(M)
mulreduce(M::Mul{<:Any, <:InfToeplitzLayouts}) = ApplyArray(M)
mulreduce(M::Mul{<:AbstractQLayout, <:InfToeplitzLayouts}) = ApplyArray(M)
simplifiable(::Mul{<:DiagonalLayout, <:InfToeplitzLayouts}) = Val(true)
simplifiable(::Mul{<:InfToeplitzLayouts, <:DiagonalLayout}) = Val(true)
mulreduce(M::Mul{<:DiagonalLayout, <:InfToeplitzLayouts}) = Lmul(M)
mulreduce(M::Mul{<:InfToeplitzLayouts, <:DiagonalLayout}) = Rmul(M)



###
# Inf-Toeplitz layout
# this could possibly be avoided via an InfFillLayout
###

const InfFill = AbstractFill{<:Any,1,<:Tuple{OneToInf}}

for Typ in (:(Tridiagonal{<:Any,<:InfFill}),
            :(SymTridiagonal{<:Any,<:InfFill}))
    @eval begin
        MemoryLayout(::Type{<:$Typ}) = TridiagonalToeplitzLayout()
        BroadcastStyle(::Type{<:$Typ}) = LazyArrayStyle{2}()
    end
end

MemoryLayout(::Type{<:Bidiagonal{<:Any,<:InfFill}}) = BidiagonalToeplitzLayout()
BroadcastStyle(::Type{<:Bidiagonal{<:Any,<:InfFill}}) = LazyArrayStyle{2}()

*(A::Bidiagonal{<:Any,<:InfFill}, B::Bidiagonal{<:Any,<:InfFill}) =
    mul(A, B)

# fall back for Ldiv
triangularlayout(::Type{<:TriangularLayout{UPLO,'N'}}, ::TridiagonalToeplitzLayout) where UPLO = BidiagonalToeplitzLayout()
materialize!(L::MatLdivVec{BidiagonalToeplitzLayout,Lay}) where Lay = materialize!(Ldiv{BidiagonalLayout{FillLayout,FillLayout},Lay}(L.A, L.B))
copyto!(dest::AbstractArray, L::Ldiv{BidiagonalToeplitzLayout,Lay}) where Lay = copyto!(dest, Ldiv{BidiagonalLayout{FillLayout,FillLayout},Lay}(L.A, L.B))


# copy for AdjOrTrans
copy(A::Adjoint{T,<:BandedMatrix{T,<:Any,OneToInf{Int}}}) where T = copy(parent(A))'
copy(A::Transpose{T,<:BandedMatrix{T,<:Any,OneToInf{Int}}}) where T = transpose(copy(parent(A)))


##
# hcat
##

Base.typed_hcat(::Type{T}, A::BandedMatrix{<:Any,<:Any,OneToInf{Int}}, B::AbstractVecOrMat...) where T = Hcat{T}(A, B...)



###
# SymTriPertToeplitz
###

MemoryLayout(::Type{<:SymTriPertToeplitz}) = PertTridiagonalToeplitzLayout()


sublayout(::ApplyBandedLayout, ::Type{<:Tuple{KR,Integer}}) where {KR<:InfAxes} =
    sublayout(PaddedColumns{UnknownLayout}(), Tuple{KR})

end # module