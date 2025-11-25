const ConstRowMatrix{T} = ApplyMatrix{T,typeof(*),<:Tuple{<:AbstractVector,<:AbstractFillMatrix{<:Any,Tuple{OneTo{Int},OneToInf{Int}}}}}
const PertConstRowMatrix{T} = Hcat{T,<:Tuple{Array{T},<:ConstRowMatrix{T}}}

struct ConstRows <: AbstractLazyLayout end
struct PertConstRows <: AbstractLazyLayout end
MemoryLayout(::Type{<:ConstRowMatrix}) = ConstRows()
MemoryLayout(::Type{<:PertConstRowMatrix}) = PertConstRows()


const TriToeplitz{T} = Tridiagonal{T,Fill{T,1,Tuple{OneToInf{Int}}}}

const SymTriPertToeplitz{T} = SymTridiagonal{T,Vcat{T,1,Tuple{Vector{T},Fill{T,1,Tuple{OneToInf{Int}}}}}}
const TriPertToeplitz{T} = Tridiagonal{T,Vcat{T,1,Tuple{Vector{T},Fill{T,1,Tuple{OneToInf{Int}}}}}}
const AdjTriPertToeplitz{T} = Adjoint{T,Tridiagonal{T,Vcat{T,1,Tuple{Vector{T},Fill{T,1,Tuple{OneToInf{Int}}}}}}}


"""
    TridiagonalToeplitzLayout

represents a matrix which is tridiagonal and toeplitz. Must support
`subdiagonalconstant`, `diagonalconstant`, `supdiagonalconstant`.
"""
struct TridiagonalToeplitzLayout <: AbstractLazyBandedLayout end

struct BidiagonalToeplitzLayout <: AbstractLazyBandedLayout end
struct PertBidiagonalToeplitzLayout <: AbstractLazyBandedLayout end
struct PertTridiagonalToeplitzLayout <: AbstractLazyBandedLayout end

const InfBaseToeplitzLayouts = Union{TridiagonalToeplitzLayout, BidiagonalToeplitzLayout, PertBidiagonalToeplitzLayout, PertTridiagonalToeplitzLayout}


subdiagonalconstant(A) = getindex_value(subdiagonaldata(A))
diagonalconstant(A) = getindex_value(diagonaldata(A))
supdiagonalconstant(A) = getindex_value(supdiagonaldata(A))

islazy_layout(::InfBaseToeplitzLayouts) = Val(true)

@inline sub_materialize(::ApplyBandedLayout{typeof(*)}, V, ::Tuple{InfAxes,InfAxes}) = V
@inline sub_materialize(::BroadcastBandedLayout, V, ::Tuple{InfAxes,InfAxes}) = V


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
    end
end

MemoryLayout(::Type{<:SymTriPertToeplitz}) = PertTridiagonalToeplitzLayout()

sublayout(::ApplyBandedLayout, ::Type{<:Tuple{KR,Integer}}) where {KR<:InfAxes} =
    sublayout(PaddedColumns{UnknownLayout}(), Tuple{KR})
sublayout(::ApplyBandedLayout, ::Type{<:Tuple{Integer,JR}}) where {JR<:InfAxes} =
    sublayout(PaddedColumns{UnknownLayout}(), Tuple{JR})

LazyArrays._vcat_sub_arguments(::Union{PertConstRows,ConstRows}, A, V) = LazyArrays._vcat_sub_arguments(ApplyLayout{typeof(hcat)}(), A, V)