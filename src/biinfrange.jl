"""
BiInfUnitRange()

Represent -∞:∞ with offset indexing
"""
struct BiInfUnitRange{T<:Real} <: AbstractInfUnitRange{T} end

BiInfUnitRange() = BiInfUnitRange{Int}()

AbstractArray{T}(a::BiInfUnitRange) where T<:Real = BiInfUnitRange{T}(a)
AbstractVector{T}(a::BiInfUnitRange) where T<:Real = BiInfUnitRange{T}(a)

unitrange(a::BiInfUnitRange) = a
Base.has_offset_axes(::BiInfUnitRange) = true

getindex(v::BiInfUnitRange{T}, i::Integer) where T = convert(T, i)
getindex(v::BiInfUnitRange{T}, i::RealInfinity) where T = i
axes(::BiInfUnitRange) = (BiInfUnitRange(),)
first(::BiInfUnitRange) = -∞
show(io::IO, ::BiInfUnitRange{Int}) = print(io, "BiInfUnitRange()")

getindex(r::BiInfUnitRange{T}, s::AbstractUnitRange{<:Integer}) where T = convert(AbstractVector{T}, s)