Array{T}(::UndefInitializer, ::NTuple{N, Infinity}) where {T,N} = throw(ArgumentError("Cannot create infinite Array"))
Array{T,N}(::UndefInitializer, ::NTuple{N, Infinity}) where {T,N} = throw(ArgumentError("Cannot create infinite Array"))
Array{T}(::UndefInitializer, ::Tuple{Integer, Infinity}) where {T,N} = throw(ArgumentError("Cannot create infinite Array"))
Array{T}(::UndefInitializer, ::Tuple{Infinity, Integer}) where {T,N} = throw(ArgumentError("Cannot create infinite Array"))
Matrix{T}(::UndefInitializer, ::Tuple{Integer, Infinity}) where T = throw(ArgumentError("Cannot create infinite Array"))
Matrix{T}(::UndefInitializer, ::Tuple{Infinity, Integer}) where T = throw(ArgumentError("Cannot create infinite Array"))

Array{T}(::UndefInitializer, ::Infinity) where T = throw(ArgumentError("Cannot create infinite Array"))

Array{T}(::UndefInitializer, ::Infinity, ::Infinity) where T = throw(ArgumentError("Cannot create infinite Array"))
Array{T}(::UndefInitializer, ::Infinity, ::Integer) where T = throw(ArgumentError("Cannot create infinite Array"))
Array{T}(::UndefInitializer, ::Integer, ::Infinity) where T = throw(ArgumentError("Cannot create infinite Array"))

Matrix{T}(::UndefInitializer, ::Infinity, ::Infinity) where T = throw(ArgumentError("Cannot create infinite Array"))
Matrix{T}(::UndefInitializer, ::Infinity, ::Integer) where T = throw(ArgumentError("Cannot create infinite Array"))
Matrix{T}(::UndefInitializer, ::Integer, ::Infinity) where T = throw(ArgumentError("Cannot create infinite Array"))

Vector{T}(::UndefInitializer, ::Infinity) where T = throw(ArgumentError("Cannot create infinite Array"))

similar(A::AbstractArray, ::Type{T}, axes::NTuple{N,OneToInf{Int}}) where {T,N} = cache(Zeros{T,N}(axes))

# This gets called when infinit number of columns
print_matrix_row(io::IO,
        X::AbstractVecOrMat, A::Vector,
        i::Integer, cols::AbstractVector{<:Infinity}, sep::AbstractString) = nothing


print_matrix_vdots(io::IO, vdots::AbstractString,
        A::Vector, sep::AbstractString, M::Integer, ::NotANumber) = nothing


# Avoid infinite loops on maximum
Base.mapreduce_impl(f, op, A::AbstractArray, ifirst::Integer, ::Infinity) =
    throw(ArgumentError("Cannot call mapreduce on an infinite length $(typeof(A))"))

function show_delim_array(io::IO, itr::AbstractArray, op, delim, cl,
                          delim_one, i1, ::Infinity)
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

BroadcastStyle(::Type{<:Eye{T,OneToInf{I}}}) where {T,I} = LazyArrayStyle{2}()

#####
# Diagonal
#####


BroadcastStyle(::Type{<:Diagonal{<:Any,<:AbstractInfUnitRange}}) = LazyArrayStyle{2}()


#####
# Vcat length
#####

function getindex(f::Vcat{T,1}, k::Infinity) where T
    length(f) == ∞ || throw(BoundsError(f,k))
    ∞
end