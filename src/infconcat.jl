# Lazy concatenation of AbstractVector's.
# Similar to Iterators.Flatten and some code has been reused from julia/base/iterators.jl

function _Vcat end

struct Vcat{T,N,I} <: AbstractArray{T,N}
    arrays::I

    global function _Vcat(A::I) where I<:Tuple{Vararg{AbstractVector{T}}} where T
        isempty(A) && throw(ArgumentError("Cannot concatenate empty vectors"))
        new{T,1,I}(A)
    end
    global function _Vcat(A::I) where I<:Tuple{Vararg{AbstractMatrix{T}}} where T
        isempty(A) && throw(ArgumentError("Cannot concatenate empty vectors"))
        m = size(A[1],2)
        for k=2:length(A)
            size(A[k],2) == m || throw(ArgumentError("number of columns of each array must match (got $(map(x->size(x,2), A)))"))
        end
        new{T,2,I}(A)
    end
end
_Vcat(::Type{T}, A) where T = _Vcat(AbstractArray{T}.(A))
_Vcat(A) = _Vcat(mapreduce(eltype, promote_type, A), A)
Vcat(args...) = _Vcat(args)
size(f::Vcat{<:Any,1}) = tuple(+(length.(f.arrays)...))
size(f::Vcat{<:Any,2}) = (+(size.(f.arrays,1)...), size(f.arrays[1],2))
Base.IndexStyle(::Type{Vcat{T,1}}) where T = Base.IndexLinear()
Base.IndexStyle(::Type{Vcat{T,2}}) where T = Base.IndexCartesian()

function getindex(f::Vcat{<:Any,1}, k::Integer)
    for A in f.arrays
        n = length(A)
        k ≤ n && return A[k]
        k -= n
    end
    throw(BoundsError("attempt to access $length(f) Vcat array."))
end

function getindex(f::Vcat{<:Any,2}, k::Integer, j::Integer)
    for A in f.arrays
        n = size(A,1)
        k ≤ n && return A[k,j]
        k -= n
    end
    throw(BoundsError("attempt to access $length(f) Vcat array."))
end

reverse(f::Vcat{<:Any,1}) = Vcat((reverse(itr) for itr in reverse(f.arrays))...)
