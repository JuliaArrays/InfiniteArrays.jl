


# This gets called when infinit number of columns
print_matrix_row(io::IO,
        X::AbstractVecOrMat, A::Vector,
        i::Integer, cols::AbstractVector{<:Infinity}, sep::AbstractString) = nothing


print_matrix_vdots(io::IO, vdots::AbstractString,
        A::Vector, sep::AbstractString, M::Integer, ::NotANumber) = nothing


# Avoid infinite loops on maximum
Base.mapreduce_impl(f, op, A::AbstractArray, ifirst::Integer, ::Infinity) =
    throw(ArgumentError("Cannot call mapreduce on an infinite length $(typeof(A))"))

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
