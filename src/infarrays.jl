


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

for typ in (:Fill, :Zeros, :Ones)
    @eval begin
        Base.IndexStyle(::Type{<:$typ{<:Any,2,Tuple{Infinity,Infinity}}}) = Base.IndexCartesian()
        Base.IndexStyle(::Type{<:$typ{<:Any,2,Tuple{Infinity,Int}}}) = Base.IndexCartesian()
    end
end
# axes(::Fill{<:Any,1,Tuple{Infinity}}) = tuple(OneToInf())

# Lazy Broadacasting
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

#####
# Diagonal
#####


BroadcastStyle(::Type{<:Diagonal{<:Any,<:AbstractInfUnitRange}}) = LazyArrayStyle{2}()


######
# PaddedArrays
######

# this is a special override that may be generalisable
broadcasted(::LazyArrayStyle{1}, op, A::Vcat{<:Any, 1, <:Tuple{<:Number, <:AbstractFill}},
                                     B::Vcat{<:Any, 1, <:Tuple{<:Number, <:AbstractFill}}) =
     Vcat(op(A.arrays[1], B.arrays[1]), op.(A.arrays[2], B.arrays[2]))

broadcasted(::LazyArrayStyle{1}, op, A::Vcat{<:Any, 1, <:Tuple{<:SVector{M}, <:AbstractFill}},
                                     B::Vcat{<:Any, 1, <:Tuple{<:SVector{M}, <:AbstractFill}}) where M =
  Vcat(op.(A.arrays[1], B.arrays[1]), op.(A.arrays[2], B.arrays[2]))
