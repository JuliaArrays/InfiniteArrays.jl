# Lazy concatenation of AbstractVector's.
# Similar to Iterators.Flatten and some code has been reused from julia/base/iterators.jl

function _Vcat end

struct Vcat{T,I} <: AbstractVector{T}
    it::I
    cumulsizes::Vector{Int}
    global _Vcat(it::I) where I<:Tuple{Vararg{AbstractVector{T}}} where T
        new{T,I}(it,cumsum(length.(it[1:end-1])))
    end
end

Vcat{T}(args...) = _Vcat(args)
length(f::Vcat) = mapreduce(+, length, f.it)

function start(f::Vcat)
    local inner, s2
    s = start(f.it)
    d = done(f.it, s)
    # this is a simple way to make this function type stable
    d && throw(ArgumentError("argument to Vcat must contain at least one iterator"))
    while !d
        inner, s = next(f.it, s)
        s2 = start(inner)
        !done(inner, s2) && break
        d = done(f.it, s)
    end
    return s, inner, s2
end

@propagate_inbounds function next(f::Vcat, state)
    s, inner, s2 = state
    val, s2 = next(inner, s2)
    while done(inner, s2) && !done(f.it, s)
        inner, s = next(f.it, s)
        s2 = start(inner)
    end
    return val, (s, inner, s2)
end

@inline function done(f::Vcat, state)
    s, inner, s2 = state
    return done(f.it, s) && done(inner, s2)
end


reverse(f::Vcat) = Vcat(reverse(itr) for itr in reverse(f.it))
