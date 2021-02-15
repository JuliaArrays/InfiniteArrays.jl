Array{T}(::UndefInitializer, ::Tuple{PosInfinity}) where T = throw(ArgumentError("Cannot create infinite Array"))
Array{T}(::UndefInitializer, ::Tuple{Integer, PosInfinity}) where {T,N} = throw(ArgumentError("Cannot create infinite Array"))
Array{T}(::UndefInitializer, ::Tuple{PosInfinity, Integer}) where {T,N} = throw(ArgumentError("Cannot create infinite Array"))
Matrix{T}(::UndefInitializer, ::Tuple{Integer, PosInfinity}) where T = throw(ArgumentError("Cannot create infinite Array"))
Matrix{T}(::UndefInitializer, ::Tuple{PosInfinity, Integer}) where T = throw(ArgumentError("Cannot create infinite Array"))
Array{T}(::UndefInitializer, ::Tuple{PosInfinity, PosInfinity}) where {T,N} = throw(ArgumentError("Cannot create infinite Array"))
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

similar(arr::AbstractArray{T}, sz::Vararg{Union{Infinity,Integer}}) where T = similar(arr, Base.to_shape(sz)...)
similar(arr::AbstractArray, ::Type{T}, sz::Vararg{Union{Infinity,Integer}}) where T = similar(arr, T, Base.to_shape(sz)...)

for Typ in (:Zeros, :Ones)
    @eval @inline $Typ{T, N}(sz::Tuple{Vararg{<:Union{Infinity,Integer}, N}}) where {T, N} = $Typ{T,N}(oneto.(sz))
end

Fill(x, n::Vararg{Union{Infinity,Integer}}) = Fill(x, Base.to_shape(n))


zeros(::Type{T}, ::Tuple{PosInfinity}) where T = cache(Zeros{T}(∞))
zeros(::Type{T}, nm::Tuple{Integer, PosInfinity}) where T = cache(Zeros{T}(nm...))
zeros(::Type{T}, nm::Tuple{PosInfinity, Integer}) where T = cache(Zeros{T}(nm...))
zeros(::Type{T}, nm::Tuple{PosInfinity, PosInfinity}) where T = cache(Zeros{T}(nm...))

ones(::Type{T}, ::Tuple{PosInfinity}) where T = cache(Ones{T}(∞))
ones(::Type{T}, nm::Tuple{Integer, PosInfinity}) where T = cache(Ones{T}(nm...))
ones(::Type{T}, nm::Tuple{PosInfinity, Integer}) where T = cache(Ones{T}(nm...))
ones(::Type{T}, nm::Tuple{PosInfinity, PosInfinity}) where T = cache(Ones{T}(nm...))

fill(x, ::Tuple{PosInfinity}) = cache(Fill(x,∞))
fill(x, nm::Tuple{Integer, PosInfinity}) = cache(Fill(x,nm...))
fill(x, nm::Tuple{PosInfinity, Integer}) = cache(Fill(x,nm...))
fill(x, nm::Tuple{PosInfinity, PosInfinity}) = cache(Fill(x,nm...))



# This gets called when infinit number of columns
axes_print_matrix_row(_, io, X, A, i, ::AbstractVector{<:PosInfinity}, sep) = nothing
print_matrix_row(io::IO, X::AbstractVecOrMat, A::Vector, i::Integer, cols::AbstractVector{<:PosInfinity}, sep::AbstractString) = nothing
Base.print_matrix_row(io::IO,
        X::Union{LayoutMatrix,
        LayoutVector,
        AbstractTriangular{<:Any,<:LayoutMatrix},
        AdjOrTrans{<:Any,<:LayoutMatrix},
        AdjOrTrans{<:Any,<:LayoutVector},
        HermOrSym{<:Any,<:LayoutMatrix},
        SubArray{<:Any,2,<:LayoutMatrix},
        Diagonal{<:Any,<:LayoutVector}}, A::Vector,
        i::Integer, cols::AbstractVector{<:PosInfinity}, sep::AbstractString) =
    axes_print_matrix_row(axes(X), io, X, A, i, cols, sep)
Base.print_matrix_row(io::IO,
        X::Union{AbstractFill{<:Any,1},
                 AbstractFill{<:Any,2},
                 Diagonal{<:Any,<:AbstractFill{<:Any,1}},
                 RectDiagonal,
                 AbstractTriangular{<:Any,<:AbstractFill{<:Any,2}}
                 }, A::Vector,
        i::Integer, cols::AbstractVector{<:PosInfinity}, sep::AbstractString) =
    axes_print_matrix_row(axes(X), io, X, A, i, cols, sep)


print_matrix_vdots(io::IO, vdots::AbstractString,
        A::Vector, sep::AbstractString, M::Integer, ::NotANumber, pad_right::Bool = true) = nothing

# Avoid infinite loops on maximum
Base.mapreduce_impl(f, op, A::AbstractArray, ifirst::Integer, ::PosInfinity) =
    throw(ArgumentError("Cannot call mapreduce on an infinite length $(typeof(A))"))

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


# Lazy Broadacasting
for typ in (:Ones, :Zeros, :Fill)
    @eval begin
        BroadcastStyle(::Type{$typ{T,N,NTuple{N,<:OneToInf}}}) where {T,N} = LazyArrayStyle{N}()
        BroadcastStyle(::Type{$typ{T,2,<:Tuple{<:Any,<:OneToInf}}}) where {T} = LazyArrayStyle{2}()
        BroadcastStyle(::Type{$typ{T,2,<:Tuple{<:OneToInf,<:Any}}}) where {T} = LazyArrayStyle{2}()
    end
end

BroadcastStyle(::Type{<:Diagonal{T,<:AbstractFill{T,1,Tuple{OneToInf{I}}}}}) where {T,I} = LazyArrayStyle{2}()

## Support broadcast(*, ::AbstractFill, A)


function copy(bc::Broadcasted{<:BroadcastStyle,<:Any,typeof(*),<:Tuple{Ones{T,1,Tuple{OneToInf{Int}}},AbstractArray{V,N}}}) where {N,T,V}
    a,b = bc.args
    @assert bc.axes == axes(b)
    convert(AbstractArray{promote_type(T,V),N}, b)
end

function copy(bc::Broadcasted{<:BroadcastStyle,<:Any,typeof(*),<:Tuple{AbstractArray{T,N},Ones{V,1,Tuple{OneToInf{Int}}}}}) where {N,T,V}
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


BroadcastStyle(::Type{<:Diagonal{<:Any,<:AbstractInfUnitRange}}) = LazyArrayStyle{2}()


#####
# Vcat length
#####

function getindex(f::Vcat{T,1}, k::PosInfinity) where T
    length(f) == ∞ || throw(BoundsError(f,k))
    ∞
end

_gettail(k, a::Number, b...) = k ≤ 1 ? tuple(a, b...) : _gettail(k - length(a), b...)
_gettail(k, a, b...) = k ≤ length(a) ? tuple(a[k:end], b...) : _gettail(k - length(a), b...)
_vcat(a) = a
_vcat(a, b, c...) = Vcat(a, b, c...)
getindex(A::Vcat, r::InfUnitRange) = Base.invoke(getindex, Tuple{AbstractArray, Any}, A, r)
_unsafe_getindex(::IndexLinear, A::Vcat, r::InfUnitRange) = _vcat(_gettail(first(r), A.args...)...)
