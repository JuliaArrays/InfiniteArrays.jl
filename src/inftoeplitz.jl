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

