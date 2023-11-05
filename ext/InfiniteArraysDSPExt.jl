module InfiniteArraysDSPExt

using InfiniteArrays
using InfiniteArrays: InfRanges, OneToInf
# Specifying the full namespace is necessary because of https://github.com/JuliaLang/julia/issues/48533
# See https://github.com/JuliaStats/LogExpFunctions.jl/pull/63
using InfiniteArrays.FillArrays
using InfiniteArrays.FillArrays: AbstractFill, getindex_value
using InfiniteArrays.LazyArrays
import DSP: conv


##
# conv
# This is useful for determining polynomial dimensions
##

conv(::Ones{T,1,<:Tuple{<:OneToInf}}, ::Ones{V,1,<:Tuple{<:OneToInf}}) where {T<:Integer,V<:Integer} =
    OneToInf{promote_type(T,V)}()
conv(::Ones{Bool,1,<:Tuple{<:OneToInf}}, ::Ones{Bool,1,<:Tuple{<:OneToInf}}) =
    OneToInf()
conv(::Ones{T,1,<:Tuple{<:OneToInf}}, ::Ones{V,1,<:Tuple{<:OneToInf}}) where {T,V} =
    one(promote_type(T,V)):∞

function conv(::Ones{T,1,<:Tuple{<:OneToInf}}, a::AbstractVector{V}) where {T,V}
    cs = cumsum(convert(AbstractVector{promote_type(T,V)}, a))
    Vcat(cs, Fill(last(cs), ∞))
end

function conv(::Ones{T,1,<:Tuple{<:OneToInf}}, a::Vector{V}) where {T,V}
    cs = cumsum(convert(AbstractVector{promote_type(T,V)}, a))
    Vcat(cs, Fill(last(cs), ∞))
end

function conv(a::AbstractVector{V}, ::Ones{T,1,<:Tuple{<:OneToInf}}) where {T,V}
    cs = cumsum(convert(AbstractVector{promote_type(T,V)}, a))
    Vcat(cs, Fill(last(cs), ∞))
end

function conv(a::Vector{V}, ::Ones{T,1,<:Tuple{<:OneToInf}}) where {T,V}
    cs = cumsum(convert(AbstractVector{promote_type(T,V)}, a))
    Vcat(cs, Fill(last(cs), ∞))
end


function conv(r::InfRanges, x::AbstractVector)
    length(x) ≠ 1 && throw(ArgumentError("conv(::$(typeof(r)), ::$(typeof(x))) not implemented"))
    first(x)*r
end
function conv(x::AbstractVector, r::InfRanges)
    length(x) ≠ 1 && throw(ArgumentError("conv(::$(typeof(r)), ::$(typeof(x))) not implemented"))
    first(x)*r
end

conv(r1::InfRanges, r2::AbstractFill{<:Any,1,<:Tuple{<:OneToInf}}) =
    cumsum(r1*getindex_value(r2))
conv(r2::AbstractFill{<:Any,1,<:Tuple{<:OneToInf}}, r1::InfRanges) =
    cumsum(getindex_value(r2)*r1)

conv(r1::InfRanges, r2::Ones{<:Any,1,<:Tuple{<:OneToInf}}) = cumsum(r1)
conv(r2::Ones{<:Any,1,<:Tuple{<:OneToInf}}, r1::InfRanges) = cumsum(r1)

conv(r1::InfRanges, r2::InfRanges) = throw(ArgumentError("conv(::$(typeof(r1)), ::$(typeof(r2))) not implemented"))

function conv(r1::AbstractFill{<:Any,1,<:Tuple{<:OneToInf}}, r2::AbstractFill{<:Any,1,<:Tuple{<:OneToInf}})
    a = getindex_value(r1) * getindex_value(r2)
    a:a:∞
end
function conv(r1::AbstractFill{<:Any,1,<:Tuple{<:OneToInf}}, r2::Ones{<:Any,1,<:Tuple{<:OneToInf}})
    a = getindex_value(r1) * getindex_value(r2)
    a:a:∞
end
function conv(r1::Ones{<:Any,1,<:Tuple{<:OneToInf}}, r2::AbstractFill{<:Any,1,<:Tuple{<:OneToInf}})
    a = getindex_value(r1) * getindex_value(r2)
    a:a:∞
end


end
