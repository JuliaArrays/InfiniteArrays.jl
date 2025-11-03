Array{T}(::UndefInitializer, ::Tuple{PosInfinity}) where T = throw(ArgumentError("Cannot create infinite Array"))
Array{T}(::UndefInitializer, ::Tuple{Integer, PosInfinity}) where T = throw(ArgumentError("Cannot create infinite Array"))
Array{T}(::UndefInitializer, ::Tuple{PosInfinity, Integer}) where T = throw(ArgumentError("Cannot create infinite Array"))
Matrix{T}(::UndefInitializer, ::Tuple{Integer, PosInfinity}) where T = throw(ArgumentError("Cannot create infinite Array"))
Matrix{T}(::UndefInitializer, ::Tuple{PosInfinity, Integer}) where T = throw(ArgumentError("Cannot create infinite Array"))
Array{T}(::UndefInitializer, ::Tuple{PosInfinity, PosInfinity}) where T = throw(ArgumentError("Cannot create infinite Array"))
Matrix{T}(::UndefInitializer, ::Tuple{PosInfinity, PosInfinity}) where T = throw(ArgumentError("Cannot create infinite Array"))

Array{T}(::UndefInitializer, ::PosInfinity) where T = throw(ArgumentError("Cannot create infinite Array"))

Array{T}(::UndefInitializer, ::PosInfinity, ::PosInfinity) where T = throw(ArgumentError("Cannot create infinite Array"))
Array{T}(::UndefInitializer, ::PosInfinity, ::Integer) where T = throw(ArgumentError("Cannot create infinite Array"))
Array{T}(::UndefInitializer, ::Integer, ::PosInfinity) where T = throw(ArgumentError("Cannot create infinite Array"))

Matrix{T}(::UndefInitializer, ::PosInfinity, ::PosInfinity) where T = throw(ArgumentError("Cannot create infinite Array"))
Matrix{T}(::UndefInitializer, ::PosInfinity, ::Integer) where T = throw(ArgumentError("Cannot create infinite Array"))
Matrix{T}(::UndefInitializer, ::Integer, ::PosInfinity) where T = throw(ArgumentError("Cannot create infinite Array"))

Vector{T}(::UndefInitializer, ::Tuple{PosInfinity}) where T = throw(ArgumentError("Cannot create infinite Array"))
Vector{T}(::UndefInitializer, ::PosInfinity) where T = throw(ArgumentError("Cannot create infinite Array"))

similar(A::AbstractArray, ::Type{T}, axes::Tuple{OneToInf{Int}}) where T = cache(Zeros{T}(axes))
similar(A::AbstractArray, ::Type{T}, axes::Tuple{OneToInf{Int},OneToInf{Int}}) where T = cache(Zeros{T}(axes))
similar(A::AbstractArray, ::Type{T}, dims::Tuple{PosInfinity}) where T = cache(Zeros{T}(dims))
similar(A::AbstractArray, ::Type{T}, dims::Tuple{PosInfinity,PosInfinity}) where T = cache(Zeros{T}(dims))
similar(A::AbstractArray, ::Type{T}, dims::Tuple{Integer,PosInfinity}) where T = cache(Zeros{T}(dims))
similar(A::AbstractArray, ::Type{T}, dims::Tuple{PosInfinity,Integer}) where T = cache(Zeros{T}(dims))

similar(::Type{<:AbstractArray{T}}, axes::Tuple{OneToInf{Int}}) where T = cache(Zeros{T}(axes))
similar(::Type{<:AbstractArray{T}}, axes::Tuple{OneToInf{Int},OneToInf{Int}}) where T = cache(Zeros{T}(axes))
similar(::Type{<:AbstractArray{T}}, axes::Tuple{OneToInf{Int},OneTo{Int}}) where T = cache(Zeros{T}(axes))
similar(::Type{<:AbstractArray{T}}, axes::Tuple{OneTo{Int},OneToInf{Int}}) where T = cache(Zeros{T}(axes))
similar(::Type{<:AbstractArray{T}}, dims::Tuple{PosInfinity}) where T = cache(Zeros{T}(dims))
similar(::Type{<:AbstractArray{T}}, dims::Tuple{PosInfinity,PosInfinity}) where T = cache(Zeros{T}(dims))
similar(::Type{<:AbstractArray{T}}, dims::Tuple{Integer,PosInfinity}) where T = cache(Zeros{T}(dims))
similar(::Type{<:AbstractArray{T}}, dims::Tuple{PosInfinity,Integer}) where T = cache(Zeros{T}(dims))

similar(arr::AbstractArray{T}, sz::Union{Infinity,Integer,AbstractUnitRange}...) where T = similar(arr, Base.to_shape(sz)...)
similar(arr::AbstractArray, ::Type{T}, sz::Union{Infinity,Integer,AbstractUnitRange}...) where T = similar(arr, T, Base.to_shape(sz)...)

for Typ in (:Zeros, :Ones)
    @eval begin
        @inline $Typ{T, N}(sz::Tuple{Vararg{Union{Infinity,Integer}, N}}) where {T, N} = $Typ{T,N}(oneto.(sz))
        @inline $Typ{T, N}(sz::Vararg{Union{Infinity,Integer}, N}) where {T, N} = $Typ{T,N}(sz)
        @inline $Typ{T}(sz::Vararg{Union{Infinity,Integer},N}) where {T, N} = $Typ{T, N}(sz)
        @inline $Typ{T}(sz::Infinity) where T = $Typ{T,1}(sz)
    end
end

Eye{T}(::Infinity) where T = Eye{T}(ℵ₀)
Eye{T}(::Infinity, ::Infinity) where T = Eye{T}(ℵ₀, ℵ₀)
Eye{T}(::Infinity, n::Integer) where T = Eye{T}(ℵ₀, n)
Eye{T}(m::Integer, ::Infinity) where T = Eye{T}(m, ℵ₀)
Eye(::Infinity) = Eye(ℵ₀)
Eye(::Infinity, ::Infinity) = Eye(ℵ₀, ℵ₀)
Eye(::Infinity, n::Integer) = Eye(ℵ₀, n)
Eye(m::Integer, ::Infinity) = Eye(m, ℵ₀)

@inline Fill(x, sz::Union{Infinity,Integer}...) = Fill(x, Base.to_shape(sz)...)
@inline Fill{T}(x, sz::Union{Infinity,Integer}...) where T = Fill{T}(x, Base.to_shape(sz)...)

zeros(sz::Union{Infinity, Integer, AbstractUnitRange}...) = zeros(Base.to_shape(sz)...)
zeros(::Type{T}, sz::Union{Infinity, Integer, AbstractUnitRange}...) where T = zeros(T, Base.to_shape(sz)...)
zeros(::Type{T}, ::Tuple{PosInfinity}) where T = cache(Zeros{T}(∞))
zeros(::Type{T}, nm::Tuple{Integer, PosInfinity}) where T = cache(Zeros{T}(nm...))
zeros(::Type{T}, nm::Tuple{Integer, Integer, PosInfinity}) where T = cache(Zeros{T}(nm...))
zeros(::Type{T}, nm::Tuple{PosInfinity, Integer}) where T = cache(Zeros{T}(nm...))
zeros(::Type{T}, nm::Tuple{PosInfinity, PosInfinity}) where T = cache(Zeros{T}(nm...))

ones(sz::Union{Infinity, Integer, AbstractUnitRange}...) = ones(Base.to_shape(sz)...)
ones(::Type{T}, ::Tuple{PosInfinity}) where T = cache(Ones{T}(∞))
ones(::Type{T}, nm::Tuple{Integer, PosInfinity}) where T = cache(Ones{T}(nm...))
ones(::Type{T}, nm::Tuple{PosInfinity, Integer}) where T = cache(Ones{T}(nm...))
ones(::Type{T}, nm::Tuple{PosInfinity, PosInfinity}) where T = cache(Ones{T}(nm...))

fill(x, sz::Union{Infinity, Integer, AbstractUnitRange}...) = fill(x, Base.to_shape(sz)...)
fill(x, ::Tuple{PosInfinity}) = cache(Fill(x,∞))
fill(x, nm::Tuple{Integer, PosInfinity}) = cache(Fill(x,nm...))
fill(x, nm::Tuple{PosInfinity, Integer}) = cache(Fill(x,nm...))
fill(x, nm::Tuple{PosInfinity, PosInfinity}) = cache(Fill(x,nm...))



# This gets called when infinite number of columns
axes_print_matrix_row(_, io, X, A, i, ::AbstractVector{<:PosInfinity}, sep) = nothing
print_matrix_row(io::IO, X::AbstractVecOrMat, A::Vector, i::Integer, cols::AbstractVector{<:PosInfinity}, sep::AbstractString, idxlast::Integer=last(axes(X, 2))) = nothing
print_matrix_row(io::IO, X::AbstractVecOrMat, A::Vector, i::Integer, cols::AbstractVector, sep::AbstractString, idxlast::Union{RealInfinity,Infinity}) = print_matrix_row(io, X, A, i, cols, sep, ℵ₀)
print_matrix_row(io::IO,
        X::Union{LayoutMatrix,
        LayoutVector,
        UpperOrLowerTriangular{<:Any,<:LayoutMatrix},
        AdjOrTrans{<:Any,<:LayoutMatrix},
        AdjOrTrans{<:Any,<:LayoutVector},
        HermOrSym{<:Any,<:LayoutMatrix},
        SubArray{<:Any,2,<:LayoutMatrix},
        Diagonal{<:Any,<:LayoutVector}}, A::Vector,
        i::Integer, cols::AbstractVector{<:PosInfinity}, sep::AbstractString, idxlast::Integer=last(axes(X, 2))) =
    axes_print_matrix_row(axes(X), io, X, A, i, cols, sep)
print_matrix_row(io::IO,
        X::Union{AbstractFill{<:Any,1},
                 AbstractFill{<:Any,2},
                 Diagonal{<:Any,<:AbstractFill{<:Any,1}},
                 RectDiagonal,
                 UpperOrLowerTriangular{<:Any,<:AbstractFill{<:Any,2}}
                 }, A::Vector,
        i::Integer, cols::AbstractVector{<:PosInfinity}, sep::AbstractString, idxlast::Integer=last(axes(X, 2))) =
    axes_print_matrix_row(axes(X), io, X, A, i, cols, sep)


print_matrix_vdots(io::IO, vdots::AbstractString,
        A::Vector, sep::AbstractString, M::Integer, ::NotANumber, pad_right::Bool = true) = nothing

# Avoid infinite loops on maximum
Base.mapreduce_impl(f, op, A::AbstractArray, ifirst::Integer, ::PosInfinity) =
    throw(ArgumentError("Cannot call mapreduce on an infinite length $(typeof(A))"))

# fix error in show
Base.isassigned(A::AbstractArray, i::RealInfinity) = i == ∞ ? isassigned(A, ℵ₀) : false
Base.getindex(A::AbstractArray, i::RealInfinity) = A[convert(Integer, i)]
Base.getindex(A::AbstractCachedVector, i::RealInfinity) = A[convert(Integer, i)]

# work around due to RealInfinity appearing from UnitStepRange
show_delim_array(io::IO, itr::AbstractArray, op, delim, cl, delim_one, i1, inf::RealInfinity) =
    show_delim_array(io, itr, op, delim, cl, delim_one, i1, convert(Integer, inf))
function show_delim_array(io::IO, itr::AbstractArray, op, delim, cl,
                          delim_one, i1, ::PosInfinity)
    print(io, op)
    l = 20
    if !show_circular(io, itr)
        recur_io = IOContext(io, :SHOWN_SET => itr)
        if !haskey(io, :compact)
          recur_io = IOContext(recur_io, :compact => true)
        end
        first = true
        i = i1
        if 20 >= i1
          while true
              if !isassigned(itr, i)
                  print(io, undef_ref_str)
              else
                  x = itr[i]
                  show(recur_io, x)
              end
              i += 1
              if i > l
                  print(io, delim)
                  print(io, ' ')
                  print(io, '…')
                  delim_one && first && print(io, delim)
                  break
              end
              first = false
              print(io, delim)
              print(io, ' ')
          end
        end
    end
    print(io, cl)
end


#####
# FillArrays
#####


# Lazy Broadcasting
for typ in (:Ones, :Zeros, :Fill)
    @eval begin
        BroadcastStyle(::Type{$typ{T,N,<:Tuple{OneToInf,Vararg{OneToInf}}}}) where {T,N} = LazyArrayStyle{N}()
        BroadcastStyle(::Type{$typ{T,2,<:Tuple{<:Any,<:OneToInf}}}) where {T} = LazyArrayStyle{2}()
        BroadcastStyle(::Type{$typ{T,2,<:Tuple{<:OneToInf,<:Any}}}) where {T} = LazyArrayStyle{2}()
    end
end

for M in (:Diagonal, :Bidiagonal, :Tridiagonal, :SymTridiagonal)
    @eval BroadcastStyle(::Type{<:$M{T,<:AbstractFill{T,1,Tuple{OneToInf{I}}}}}) where {T,I} = LazyArrayStyle{2}()
end

## Support broadcast(*, ::AbstractFill, A)


function copy(bc::Broadcasted{<:BroadcastStyle,<:Any,typeof(*),<:Tuple{Ones{T,1,Tuple{OneToInf{Int}}},AbstractArray{V,N}}}) where {N,T,V}
    a,b = bc.args
    @assert bc.axes == axes(b)
    convert(AbstractArray{promote_type(T,V),N}, b)
end

function copy(bc::Broadcasted{<:BroadcastStyle,<:Any,typeof(*),<:Tuple{Ones{T,2,Tuple{OneToInf{Int},OneTo{Int}}},AbstractArray{V,N}}}) where {N,T,V}
    a,b = bc.args
    @assert bc.axes == axes(b)
    convert(AbstractArray{promote_type(T,V),N}, b)
end

function _onesbroadcast_ifinf(::NTuple{2,OneToInf{Int}}, bc)
    a,b = bc.args
    convert(AbstractArray{promote_type(eltype(a),eltype(b))}, b)
end
_onesbroadcast_ifinf(_, bc) = BroadcastArray(bc)

function copy(bc::Broadcasted{<:BroadcastStyle,<:Any,typeof(*),<:Tuple{Ones{T,2,NTuple{2,OneToInf{Int}}},AbstractArray{V,N}}}) where {N,T,V}
    a,b = bc.args
    _onesbroadcast_ifinf(axes(b), bc)
end



copy(bc::Broadcasted{<:BroadcastStyle,<:Any,typeof(*),<:Tuple{Ones{T,1,Tuple{OneToInf{Int}}},AbstractArray{V,N},Vararg{Any}}}) where {N,T,V} =
    broadcast(*, first(bc.args), broadcast(*, tail(bc.args)...))
copy(bc::Broadcasted{<:BroadcastStyle,<:Any,typeof(*),<:Tuple{Ones{T,2,Tuple{OneToInf{Int},OneTo{Int}}},AbstractArray{V,N},Vararg{Any}}}) where {N,T,V} =
    broadcast(*, first(bc.args), broadcast(*, tail(bc.args)...))
copy(bc::Broadcasted{<:BroadcastStyle,<:Any,typeof(*),<:Tuple{Ones{T,2,NTuple{2,OneToInf{Int}}},AbstractArray{V,N},Vararg{Any}}}) where {N,T,V} =
    broadcast(*, first(bc.args), broadcast(*, tail(bc.args)...))    

function copy(bc::Broadcasted{<:BroadcastStyle,<:Any,typeof(*),<:Tuple{AbstractArray{T,N},Ones{V,1,Tuple{OneToInf{Int}}}}}) where {N,T,V}
    a,b = bc.args
    @assert bc.axes == axes(a)
    convert(AbstractArray{promote_type(T,V),N}, a)
end

function copy(bc::Broadcasted{<:BroadcastStyle,<:Any,typeof(*),<:Tuple{AbstractArray{T,N},Ones{V,2,Tuple{OneTo{Int},OneToInf{Int}}}}}) where {N,T,V}
    a,b = bc.args
    @assert bc.axes == axes(a)
    convert(AbstractArray{promote_type(T,V),N}, a)
end

function copy(bc::Broadcasted{<:BroadcastStyle,<:Any,typeof(*),<:Tuple{AbstractArray{T,N},Ones{V,2,Tuple{OneToInf{Int},OneToInf{Int}}}}}) where {N,T,V}
    a,b = bc.args
    @assert bc.axes == axes(a)
    convert(AbstractArray{promote_type(T,V),N}, a)
end

function copy(bc::Broadcasted{<:BroadcastStyle,<:Any,typeof(*),<:Tuple{AbstractFill{T,1,Tuple{OneToInf{Int}}},AbstractArray{V,N}}}) where {N,T<:Number,V}
    a,b = bc.args
    @assert bc.axes == axes(b)
    getindex_value(a) * b
end

function copy(bc::Broadcasted{<:BroadcastStyle,<:Any,typeof(*),<:Tuple{AbstractArray{T,N},AbstractFill{V,1,Tuple{OneToInf{Int}}}}}) where {N,T,V<:Number}
    a,b = bc.args
    @assert bc.axes == axes(a)
    a * getindex_value(b)
end

function copy(bc::Broadcasted{<:BroadcastStyle,<:Any,typeof(*),<:Tuple{AbstractFill{T,1,Tuple{OneToInf{Int}}},AbstractArray{V,N}}}) where {N,T,V}
    a,b = bc.args
    @assert bc.axes == axes(b)
    Ref(getindex_value(a)) .* b # Use broadcast in-case a is array-valued
end

function copy(bc::Broadcasted{<:BroadcastStyle,<:Any,typeof(*),<:Tuple{AbstractArray{T,N},AbstractFill{V,1,Tuple{OneToInf{Int}}}}}) where {N,T,V}
    a,b = bc.args
    @assert bc.axes == axes(a)
    a .* Ref(getindex_value(b)) # Use broadcast in-case b is array-valued
end

# row Vector
function copy(bc::Broadcasted{<:BroadcastStyle,<:Any,typeof(*),<:Tuple{Adjoint{T,Ones{T,1,Tuple{OneToInf{Int}}}},AbstractMatrix{V}}}) where {T,V}
    a,b = bc.args
    @assert bc.axes == axes(b)
    convert(AbstractMatrix{promote_type(T,V)}, b)
end

function copy(bc::Broadcasted{<:BroadcastStyle,<:Any,typeof(*),<:Tuple{AbstractMatrix{T},Adjoint{V,Ones{V,1,Tuple{OneToInf{Int}}}}}}) where {T,V}
    a,b = bc.args
    @assert bc.axes == axes(a)
    convert(AbstractMatrix{promote_type(T,V)}, a)
end

function copy(bc::Broadcasted{<:BroadcastStyle,<:Any,typeof(*),<:Tuple{Adjoint{T,<:AbstractFill{T,1,Tuple{OneToInf{Int}}}},AbstractMatrix{V}}}) where {T,V}
    a,b = bc.args
    @assert bc.axes == axes(b)
    getindex_value(a') * b
end

function copy(bc::Broadcasted{<:BroadcastStyle,<:Any,typeof(*),<:Tuple{AbstractMatrix{T},Adjoint{V,<:AbstractFill{V,1,Tuple{OneToInf{Int}}}}}}) where {T,V}
    a,b = bc.args
    @assert bc.axes == axes(a)
    a * getindex_value(b')
end


#####
# Diagonal
#####

one(D::Diagonal{T,<:AbstractFill{T,1,Tuple{OneToInf{Int}}}}) where T = Eye{T}(size(D,1))
copy(D::Diagonal{T,<:AbstractFill{T,1,Tuple{OneToInf{Int}}}}) where T = D

for M in (:Diagonal, :Bidiagonal, :Tridiagonal, :SymTridiagonal)
    @eval BroadcastStyle(::Type{<:$M{<:Any,<:AbstractInfUnitRange}}) = LazyArrayStyle{2}()
end
sub_materialize(::AbstractBandedLayout, V, ::Tuple{InfAxes,InfAxes}) = V
sub_materialize(::AbstractBandedLayout, V, ::Tuple{OneTo{Int},InfAxes}) = V
sub_materialize(::AbstractBandedLayout, V, ::Tuple{InfAxes,OneTo{Int}}) = V

##
# banded columns are padded
##

sublayout(::DiagonalLayout{L}, ::Type{<:Tuple{KR,Integer}}) where {L,KR<:InfAxes} =
    sublayout(PaddedColumns{UnknownLayout}(), Tuple{KR})
sublayout(::DiagonalLayout{L}, ::Type{<:Tuple{Integer,JR}}) where {L,JR<:InfAxes} =
    sublayout(PaddedColumns{UnknownLayout}(), Tuple{JR})
-

sublayout(::AbstractBandedLayout, ::Type{<:Tuple{KR,Integer}}) where {KR<:InfAxes} =
    sublayout(PaddedColumns{UnknownLayout}(), Tuple{KR})
sublayout(::AbstractBandedLayout, ::Type{<:Tuple{Integer,JR}}) where {JR<:InfAxes} =
    sublayout(PaddedColumns{UnknownLayout}(), Tuple{JR})


function sub_paddeddata(::AbstractBandedLayout, S::SubArray{T,1,<:AbstractMatrix,<:Tuple{InfAxes,Integer}}) where T
    P = parent(S)
    (kr,j) = parentindices(S)
    view(P,first(kr):last(colsupport(P,j)),j)
end

function sub_paddeddata(::AbstractBandedLayout, S::SubArray{T,1,<:AbstractMatrix,<:Tuple{Integer,InfAxes}}) where T
    P = parent(S)
    (k,jr) = parentindices(S)
    view(P,k,first(jr):last(rowsupport(P,k)))
end


#####
# Concat: Hcat/Vcat
#####

function getindex(f::Vcat{T,1}, k::PosInfinity) where T
    length(f) == ℵ₀ || throw(BoundsError(f,k))
    ∞
end

_gettail(k, a::Number, b...) = k ≤ 1 ? tuple(a, b...) : _gettail(k - length(a), b...)
_gettail(k, a, b...) = k ≤ length(a) ? tuple(a[k:end], b...) : _gettail(k - length(a), b...)
_vcat(a) = a
_vcat(a, b, c...) = Vcat(a, b, c...)
getindex(A::Vcat, r::InfUnitRange) = Base.invoke(getindex, Tuple{AbstractArray, Any}, A, r)
_unsafe_getindex(::IndexLinear, A::Vcat, r::InfUnitRange) = _vcat(_gettail(first(r), A.args...)...)

# some common cases not caught by LayoutArrays + ambiguities
for InfColMatrix in (:(SubArray{<:Any,2,<:Any,<:Tuple{Any,InfIndexRanges}}),
                     :(SubArray{<:Any,2,<:LayoutVecOrMat,<:Tuple{Any,InfIndexRanges}}),
                     :(AbstractFill{<:Any,2,Tuple{OneTo{Int},OneToInf{Int}}}))
    @eval begin
        Base.typed_vcat(::Type{T}, A::$InfColMatrix, B::AbstractVecOrMat...) where T = Vcat{T}(A, B...)
        Base.typed_vcat(::Type{T}, A::$InfColMatrix, B::LayoutVecOrMats, C::AbstractVecOrMat...) where T = Vcat{T}(A, B, C...)
        Base.hcat(A::Number, B::$InfColMatrix) = Hcat(A, B)
        Base.typed_hcat(::Type{T}, A::AbstractVecOrMat, B::$InfColMatrix) where T = Hcat{T}(A, B)
    end
end

for InfRowArray in (:(AbstractFill{<:Any,1,Tuple{OneToInf{Int}}}),
                    :(AbstractFill{<:Any,2,Tuple{OneToInf{Int},OneTo{Int}}}),
                    :(SubArray{<:Any,2,<:Any,<:Tuple{InfIndexRanges,Any}}),
                    :(SubArray{<:Any,2,<:LayoutVecOrMat,<:Tuple{InfIndexRanges,Any}}),
                    :(SubArray{<:Any,1,<:Any,<:Tuple{Any,InfIndexRanges}}),
                    :(SubArray{<:Any,1,<:LayoutMatrix,<:Tuple{Any,InfIndexRanges}}),
                    :(SubArray{<:Any,1,<:Any,<:Tuple{InfIndexRanges,Any}}),
                    :(SubArray{<:Any,1,<:LayoutMatrix,<:Tuple{InfIndexRanges,Any}}))
    @eval begin
        Base.typed_hcat(::Type{T}, A::$InfRowArray, B::AbstractVecOrMat...) where T = Hcat{T}(A, B...)
        Base.typed_hcat(::Type{T}, A::$InfRowArray, B::LayoutVecOrMats, C::AbstractVecOrMat...) where T = Hcat{T}(A, B, C...)
        Base.typed_hcat(::Type{T}, A::AbstractVecOrMat, B::$InfRowArray, C::AbstractVecOrMat...) where T = Hcat{T}(A, B, C...)
        Base.typed_hcat(::Type{T}, A::LayoutVecOrMat, B::$InfRowArray, C::AbstractVecOrMat...) where T = Hcat{T}(A, B, C...)
    end
end


ArrayLayouts.typed_hcat(::Type{T}, ::Tuple{InfiniteCardinal{0},Any}, A...) where T = Hcat{T}(A...)
ArrayLayouts.typed_hcat(::Type{T}, ::Tuple{InfiniteCardinal{0},InfiniteCardinal{0}}, A...) where T = Hcat{T}(A...)
ArrayLayouts.typed_hcat(::Type{T}, ::Tuple{Any,InfiniteCardinal{0}}, A...) where T = Hcat{T}(A...)

ArrayLayouts.typed_vcat(::Type{T}, ::Tuple{InfiniteCardinal{0},Any}, A...) where T = Vcat{T}(A...)
ArrayLayouts.typed_vcat(::Type{T}, ::Tuple{InfiniteCardinal{0},InfiniteCardinal{0}}, A...) where T = Vcat{T}(A...)
ArrayLayouts.typed_vcat(::Type{T}, ::Tuple{Any,InfiniteCardinal{0}}, A...) where T = Vcat{T}(A...)

##
# lazy sub_materialize
##

sub_materialize(_, V, ::Tuple{InfAxes}) = V
sub_materialize(_, V, ::Tuple{InfAxes,InfAxes}) = V
sub_materialize(_, V, ::Tuple{Any,InfAxes}) = V
sub_materialize(_, V, ::Tuple{InfAxes,Any}) = V


sub_materialize(::ApplyLayout{typeof(vcat)}, V::AbstractVector, ::Tuple{InfAxes}) = ApplyArray(V)
sub_materialize(::ApplyLayout{typeof(vcat)}, V::AbstractMatrix, ::Tuple{InfAxes, InfAxes}) = ApplyArray(V)
sub_materialize(::ApplyLayout{typeof(vcat)}, V::AbstractMatrix, ::Tuple{Any, InfAxes}) = ApplyArray(V)
sub_materialize(::ApplyLayout{typeof(vcat)}, V::AbstractMatrix, ::Tuple{InfAxes, Any}) = ApplyArray(V)

sub_materialize(::ApplyLayout{typeof(hcat)}, V, ::Tuple{InfAxes}) = V
sub_materialize(::ApplyLayout{typeof(hcat)}, V::AbstractMatrix, ::Tuple{InfAxes, InfAxes}) = ApplyArray(V)
sub_materialize(::ApplyLayout{typeof(hcat)}, V::AbstractMatrix, ::Tuple{Any, InfAxes}) = ApplyArray(V)
sub_materialize(::ApplyLayout{typeof(hcat)}, V::AbstractMatrix, ::Tuple{InfAxes, Any}) = ApplyArray(V)


sub_materialize(::AbstractPaddedLayout, v::AbstractVector, ::Tuple{InfAxes}) = _padded_sub_materialize(v)

sub_materialize(::AbstractPaddedLayout, v::AbstractMatrix, ::Tuple{InfAxes, InfAxes}) = v
sub_materialize(::AbstractPaddedLayout, v::AbstractMatrix, ::Tuple{InfAxes, Any}) = v
sub_materialize(::AbstractPaddedLayout, v::AbstractMatrix, ::Tuple{Any, InfAxes}) = v

sub_materialize(lay::InvColumnLayout, v::AbstractVector, ax::Tuple{InfAxes}) =
    Base.invoke(sub_materialize, Tuple{InvColumnLayout, AbstractVector, Any}, lay, v, ax)


Base._unsafe_getindex(::IndexStyle, A::AbstractVector, r::InfAxes) = layout_getindex(A, r)
Base._unsafe_getindex(::IndexStyle, A::AbstractFill{<:Any,1}, r::InfAxes) = FillArrays._fill_getindex(A, r)
getindex(A::AbstractCachedVector, r::InfAxes) = layout_getindex(A, r)
# preserve padded/fill structure
getindex(A::CachedVector{<:Any,<:AbstractVector,<:AbstractFill{<:Any,1}}, r::InfAxes) = LazyArrays.cache_getindex(A, r)
# don't resize to ∞
Base.isassigned(A::AbstractCachedVector, r::InfiniteCardinal{0}) = true
getindex(A::AbstractCachedVector, r::InfiniteCardinal{0}) = A.array[r]



Base._unsafe_getindex(::IndexStyle, A::AbstractMatrix, kr::InfAxes, jr::InfAxes) = layout_getindex(A, kr, jr)
Base._unsafe_getindex(::IndexStyle, A::AbstractMatrix, kr::Union{Real, AbstractArray}, jr::InfAxes) = layout_getindex(A, kr, jr)
Base._unsafe_getindex(::IndexStyle, A::AbstractMatrix, kr::InfAxes, jr::Union{Real, AbstractArray}) = layout_getindex(A, kr, jr)

Base._unsafe_getindex(::IndexStyle, A::AbstractFill{<:Any,2}, kr::InfAxes, jr::InfAxes) = FillArrays._fill_getindex(A, kr, jr)
Base._unsafe_getindex(::IndexStyle, A::AbstractFill{<:Any,2}, kr::Union{Real, AbstractArray}, jr::InfAxes) = FillArrays._fill_getindex(A, kr, jr)
Base._unsafe_getindex(::IndexStyle, A::AbstractFill{<:Any,2}, kr::InfAxes, jr::Union{Real, AbstractArray}) = FillArrays._fill_getindex(A, kr, jr)

@inline getindex(A::ApplyMatrix{<:Any,typeof(hcat)}, kr::InfAxes, j::Integer) = layout_getindex(A, kr, j)

Base.checkindex(::Type{Bool}, inds::AbstractInfUnitRange, I::AbstractFill) = Base.checkindex(Bool, inds, getindex_value(I))
Base.checkindex(::Type{Bool}, inds::AbstractInfUnitRange, I::AbstractFill{Bool}) = axes(I,1) == inds
LazyArrays.cache_getindex(::InfiniteCardinal{0}, A::AbstractVector, I, J...) = layout_getindex(A, I, J...)
LazyArrays.cache_getindex(::InfiniteCardinal{0}, A::CachedVector{<:Any,<:AbstractVector,<:AbstractFill{<:Any,1}}, I::AbstractVector) = LazyArrays.cache_getindex(nothing, A, I)


*(a::AbstractVector, b::AbstractFill{<:Any,2,Tuple{OneTo{Int},OneToInf{Int}}}) = ApplyArray(*,a,b)
