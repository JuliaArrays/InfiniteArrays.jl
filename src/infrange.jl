# This file is mmodified from Julia. License is MIT: https://julialang.org/license

(:)(start::T, stop::Infinity) where {T<:Integer} = InfUnitRange{T}(start)
function (:)(start::T, step::T, stop::OrientedInfinity) where {T<:Real}
    sign(step) == sign(stop) || throw(ArgumentError("InfStepRange must have infinite length"))
    InfStepRange(start, step)
end
(:)(start::T, step::Real, stop::OrientedInfinity) where {T<:Real} = (:)(promote(start, step)..., stop)
(:)(start::Real, step, stop::Infinity)= (:)(start, step, OrientedInfinity(stop))
(:)(::Infinity, _, ::Real) = throw(ArgumentError("Cannot create range starting at infinity"))

# AbstractFloat specializations
(:)(a::T, b::Union{Infinity,OrientedInfinity}) where {T<:Real} = (:)(a, T(1), b)

function (:)(start::T, step::T, stop::Infinity) where {T<:Real}
    sign(step) == sign(stop) || throw(ArgumentError("InfStepRange must have infinite length"))
    InfStepRange(start,step)
end

# this is needed for showarray
(:)(::Infinity, ::Infinity) = 1:0


# Range of a given length: range(a, [step=s,] length=l), no stop
_range(a::Real,          ::Nothing,         ::Nothing, len::Infinity) = InfUnitRange{typeof(a)}(a)
_range(a::AbstractFloat, ::Nothing,         ::Nothing, len::Infinity) = _range(a, oftype(a, 1),   nothing, len)
_rangestyle(::Ordered, ::ArithmeticWraps, a::T, step::S, len::Infinity) where {T,S} =
    InfStepRange{T,S}(a, step)
_range(a::T, st::T, ::Nothing, ::Infinity) where T<:Union{Float16,Float32,Float64} =
    InfStepRange{T,T}(a, st)

# Construct range for rational start=start_n/den, step=step_n/den
floatrange(::Type{T}, start_n::Integer, step_n::Integer, ::Infinity, den::Integer) where T =
    InfStepRange(T(start_n)/den,T(step_n)/den)

floatrange(a::AbstractFloat, st::AbstractFloat, ::Infinity, divisor::AbstractFloat) =
    InfStepRange(a/divisor,st/divisor)



## 1-dimensional ranges ##



struct InfStepRange{T,S} <: OrdinalRange{T,S}
    start::T
    step::S
    function InfStepRange{T,S}(start::T, step::S) where {T,S}
        new{T,S}(start, step)
    end
end

InfStepRange(start::T, step::S) where {T,S} = InfStepRange{T,S}(start,step)
InfStepRange{T,S}(start, step) where {T,S} = InfStepRange{T,S}(convert(T,start),convert(S,step))

abstract type AbstractInfUnitRange{T<:Real} <: AbstractUnitRange{T} end

done(r::AbstractInfUnitRange{T}, i) where {T} = false
unitrange_last(start, stop::Infinity) = ∞

struct InfUnitRange{T<:Real} <: AbstractInfUnitRange{T}
    start::T
end

InfUnitRange(a::InfUnitRange) = a
InfUnitRange{T}(a::AbstractInfUnitRange) where T<:Real = InfUnitRange{T}(first(a))
AbstractArray{T}(a::InfUnitRange) where T<:Real = InfUnitRange{T}(a.start)
AbstractVector{T}(a::InfUnitRange) where T<:Real = InfUnitRange{T}(a.start)



"""
    OneToInf(n)

Define an `AbstractInfUnitRange` that behaves like `1:∞`, with the added
distinction that the limits are guaranteed (by the type system) to
be 1 and ∞.
"""
struct OneToInf{T<:Integer} <: AbstractInfUnitRange{T} end

OneToInf() = OneToInf{Int}()

AbstractArray{T}(a::OneToInf) where T<:Integer = OneToInf{T}()
AbstractVector{T}(a::OneToInf) where T<:Integer = OneToInf{T}()
AbstractArray{T}(a::OneToInf) where T<:Real = InfUnitRange{T}(a)
AbstractVector{T}(a::OneToInf) where T<:Real = InfUnitRange{T}(a)


(==)(::OneToInf, ::OneToInf) = true

## interface implementations

const InfRanges{T} = Union{InfStepRange{T},AbstractInfUnitRange{T}}

size(r::InfRanges) = (∞,)

isempty(r::InfRanges) = false

step(r::InfStepRange) = r.step

length(r::InfRanges) = ∞
unsafe_length(r::InfRanges) = ∞

first(r::OneToInf{T}) where {T} = oneunit(T)

last(r::AbstractInfUnitRange) = ∞
last(r::InfStepRange) = sign(step(r))*∞

minimum(r::InfUnitRange) = first(r)
maximum(r::InfUnitRange) = last(r)


## iteration
start(r::InfStepRange) = oftype(r.start + r.step, r.start)
next(r::InfStepRange{T}, i) where {T} = (convert(T,i), i+r.step)

start(r::InfUnitRange{T}) where {T} = oftype(r.start + oneunit(T), r.start)
start(r::OneToInf{T}) where {T} = oneunit(T)

done(r::InfStepRange{T}, i) where {T} = false

## indexing

_sub2ind(inds::Tuple{OneToInf}, i::Integer)    = i


function getindex(v::InfUnitRange{T}, i::Integer) where T
    @boundscheck i > 0 || Base.throw_boundserror(v, i)
    convert(T, first(v) + i - 1)
end
function getindex(v::OneToInf{T}, i::Integer) where T
    @boundscheck i > 0 || Base.throw_boundserror(v, i)
    convert(T, i)
end
function getindex(v::InfStepRange{T}, i::Integer) where T
    @boundscheck i > 0 || Base.throw_boundserror(v, i)
    convert(T, first(v) + (i - 1)*step(v))
end

getindex(::AbstractInfUnitRange, ::Infinity) = ∞
getindex(::OneToInf, ::Infinity) = ∞
getindex(v::InfUnitRange{T}, i::Infinity) where T = ∞
getindex(v::OneToInf{T}, i::Infinity) where T = ∞
getindex(v::InfStepRange{T}, i::Infinity) where T = ∞

function getindex(r::AbstractInfUnitRange, s::AbstractInfUnitRange{<:Integer})
    f = first(r)
    @boundscheck first(s) ≥ 1 || Base.throw_boundserror(r, first(s))
    st = oftype(f, f + first(s)-1)
    st:∞
end

function getindex(r::AbstractInfUnitRange, s::AbstractUnitRange{<:Integer})
    f = first(r)
    @boundscheck first(s) ≥ 1 || Base.throw_boundserror(r, first(s))
    st = oftype(f, f + first(s)-1)
    range(st; length=length(s))
end

getindex(r::OneToInf{T}, s::OneTo) where T = OneTo(T(s.stop))

function getindex(r::AbstractInfUnitRange, s::InfStepRange{<:Integer})
    @_inline_meta
    @boundscheck (step(s) > 0 && first(s) ≥ 1) || throw(BoundsError(r, minimum(s)))
    st = oftype(first(r), first(r) + s.start-1)
    new_step = step(s)
    st:new_step:sign(new_step)*∞
end

function getindex(r::AbstractInfUnitRange, s::StepRange{<:Integer})
    @_inline_meta
    @boundscheck minimum(s) ≥ 1 || throw(BoundsError(r, minimum(s)))
    st = oftype(first(r), first(r) + s.start-1)
    range(st; step=step(s), length=length(s))
end

function getindex(r::InfStepRange, s::InfRanges{<:Integer})
    @_inline_meta
    @boundscheck (step(s) > 0 && first(s) ≥ 1) || throw(BoundsError(r, minimum(s)))
    st = oftype(r.start, r.start + (first(s)-1)*step(r))
    new_step = step(r)*step(s)
    st:new_step:sign(new_step)*∞
end

function getindex(r::InfStepRange, s::AbstractRange{<:Integer})
    @_inline_meta
    @boundscheck isempty(s) || minimum(s) ≥ 1 || throw(BoundsError(r, minimum(s)))
    st = oftype(r.start, r.start + (first(s)-1)*step(r))
    range(st; step=step(r)*step(s), length=length(s))
end

show(io::IO, r::InfUnitRange) = print(io, repr(first(r)), ':', repr(last(r)))
show(io::IO, r::OneToInf{Int}) = print(io, "OneToInf()")
show(io::IO, r::OneToInf{T}) where T = print(io, "OneToInf{$T}()")

intersect(r::OneToInf{T}, s::OneToInf{V}) where {T,V} = OneToInf{promote_type(T,V)}()
intersect(r::OneToInf{T}, s::OneTo{T}) where T = s
intersect(r::OneTo{T}, s::OneToInf{T}) where T = r

intersect(r::AbstractInfUnitRange{<:Integer}, s::AbstractInfUnitRange{<:Integer}) =
    InfUnitRange(max(first(r),first(s)))
intersect(r::AbstractInfUnitRange{<:Integer}, s::AbstractUnitRange{<:Integer}) =
    max(first(r),first(s)):last(s)
intersect(r::AbstractUnitRange{<:Integer}, s::AbstractInfUnitRange{<:Integer}) =
    intersect(s, r)


intersect(i::Integer, r::AbstractInfUnitRange{<:Integer}) =
    i < first(r) ? (first(r):i) : (i:i)

intersect(r::AbstractInfUnitRange{<:Integer}, i::Integer) = intersect(i, r)

function intersect(r::AbstractInfUnitRange{<:Integer}, s::StepRange{<:Integer})
    if isempty(s)
        range(first(r); length=0)
    elseif step(s) == 0
        intersect(first(s), r)
    elseif step(s) < 0
        intersect(r, reverse(s))
    else
        sta = first(s)
        ste = step(s)
        sto = last(s)
        lo = first(r)
        i0 = max(sta, lo + mod(sta - lo, ste))
        i1 = sto
        i0:ste:i1
    end
end

intersect(r::StepRange{<:Integer}, s::AbstractInfUnitRange{<:Integer}) =
    intersect(s, r)


function intersect(r::InfStepRange, s::InfStepRange)
    sign(step(r)) == sign(step(s)) || throw(ArgumentError("Can only take intersection of infinite ranges"))
    start1 = first(r)
    step1 = step(r)
    start2 = first(s)
    step2 = step(s)
    a = lcm(step1, step2)

    g, x, y = gcdx(step1, step2)

    if rem(start1 - start2, g) != 0
        # Unaligned, no overlap possible.
        throw(ArgumentError("Cannot take intersection of InfStepRange that has no elements"))
    end

    z = div(start1 - start2, g)
    b = start1 - x * z * step1
    # Possible points of the intersection of r and s are
    # ..., b-2a, b-a, b, b+a, b+2a, ...
    # Determine where in the sequence to start and stop.
    m = max(start1 + mod(b - start1, a), start2 + mod(b - start2, a))
    m:a:∞
end

intersect(r::InfStepRange, s::AbstractInfUnitRange) = intersect(r, InfStepRange(s))
intersect(r::AbstractInfUnitRange, s::InfStepRange) = intersect(InfStepRange(r), s)

function intersect(r::InfStepRange, s::StepRange)
    if isempty(s)
        return range(first(r); step=step(r), length=0)
    elseif step(s) < 0
        return intersect(r, reverse(s))
    end

    start1 = first(r)
    step1 = step(r)
    start2 = first(s)
    step2 = step(s)
    stop2 = last(s)
    a = lcm(step1, step2)

    g, x, y = gcdx(step1, step2)

    if rem(start1 - start2, g) != 0
        # Unaligned, no overlap possible.
        return range(start1; step=a, length=0)
    end

    z = div(start1 - start2, g)
    b = start1 - x * z * step1
    # Possible points of the intersection of r and s are
    # ..., b-2a, b-a, b, b+a, b+2a, ...
    # Determine where in the sequence to start and stop.
    m = max(start1 + mod(b - start1, a), start2 + mod(b - start2, a))
    n = stop2 - mod(stop2 - b, a)
    m:a:n
end

intersect(s::StepRange, r::InfStepRange) = intersect(r, s)
intersect(s::AbstractRange, r::InfStepRange) = intersect(StepRange(s), r)
intersect(s::InfStepRange, r::AbstractRange) = intersect(s, StepRange(r))


promote_rule(a::Type{InfUnitRange{T1}}, b::Type{InfUnitRange{T2}}) where {T1,T2} =
    InfUnitRange{promote_type(T1,T2)}
InfUnitRange{T}(r::InfUnitRange{T}) where {T<:Real} = r
InfUnitRange{T}(r::UnitRange) where {T<:Real} = InfUnitRange{T}(r.start)

promote_rule(a::Type{OneToInf{T1}}, b::Type{OneToInf{T2}}) where {T1,T2} =
    OneToInf{promote_type(T1,T2)}
OneToInf{T}(r::OneToInf{T}) where {T<:Integer} = r
OneToInf{T}(r::OneToInf) where {T<:Integer} = OneToInf{T}()


promote_rule(::Type{InfStepRange{T1a,T1b}}, ::Type{InfStepRange{T2a,T2b}}) where {T1a,T1b,T2a,T2b} =
    InfStepRange{promote_type(T1a,T2a), promote_type(T1b,T2b)}
promote_rule(a::Type{InfStepRange{T1a,T1b}}, ::Type{UR}) where {T1a,T1b,UR<:AbstractInfUnitRange} =
    promote_rule(a, InfStepRange{eltype(UR), eltype(UR)})

InfStepRange{T1,T2}(r::InfStepRange{T1,T2}) where {T1,T2} = r

InfStepRange{T1,T2}(r::AbstractRange) where {T1,T2} =
    InfStepRange{T1,T2}(convert(T1, first(r)), convert(T2, step(r)))
InfStepRange(r::InfUnitRange{T}) where {T} =
    InfStepRange{T,T}(first(r), step(r))
(::Type{InfStepRange{T1,T2} where T1})(r::AbstractRange) where {T2} =
    InfStepRange{eltype(r),T2}(r)





## sorting ##

sum(r::InfRanges{<:Real}) = last(r)
mean(r::InfRanges{<:Real}) = last(r)
median(r::InfRanges{<:Real}) = last(r)

in(x::Union{Infinity,OrientedInfinity}, r::InfRanges) = false # never reach it...
in(x::Infinity, r::InfRanges{<:Integer}) = false # never reach it...
in(x::Real, r::InfRanges{<:Real}) = _in_range(x, r)
# This method needs to be defined separately since -(::T, ::T) can be implemented
# even if -(::T, ::Real) is not
in(x::T, r::InfRanges{T}) where {T} = _in_range(x, r)

in(x::Integer, r::AbstractInfUnitRange{<:Integer}) = (first(r) <= x)

in(x::Real, r::InfRanges{T}) where {T<:Integer} =
    isinteger(x) && !isempty(r) && ifelse(step(r) > 0, x ≥ first(r), x ≤ first(r)) &&
            (mod(convert(T,x),step(r))-mod(first(r),step(r)) == 0)

# Addition/subtraction of ranges
-(r1::OneToInf{T}, r2::OneToInf{V}) where {T,V} = Zeros{promote_type(T,V)}(∞)
-(r1::AbstractInfUnitRange, r2::AbstractInfUnitRange) = Fill(first(r1)-first(r2), ∞)




### Lazy broadcasting

BroadcastStyle(::Type{<:InfRanges}) = LazyArrayStyle{1}()



const InfIndexRanges{T<:Integer} = Union{InfStepRange{T},AbstractInfUnitRange{T},Slice{OneToInf{T}}}

BroadcastStyle(::Type{<:SubArray{<:Any,1,<:Any,Tuple{<:InfIndexRanges}}})= LazyArrayStyle{1}()
BroadcastStyle(::Type{<:SubArray{<:Any,2,<:Any,<:Tuple{<:InfIndexRanges,<:InfIndexRanges}}})= LazyArrayStyle{1}()
BroadcastStyle(::Type{<:SubArray{<:Any,2,<:Any,<:Tuple{<:InfIndexRanges,<:Any}}})= LazyArrayStyle{1}()
BroadcastStyle(::Type{<:SubArray{<:Any,2,<:Any,<:Tuple{<:Any,<:InfIndexRanges}}})= LazyArrayStyle{1}()


cumsum(r::InfRanges) = OneToInf() .* (first(r) .+ r) .÷ 2

diff(r::InfRanges) = Fill(step(r),∞)
diff(r::OneToInf{T}) where T = Ones{T}(∞)


##
# conv
# This is useful for determining polynomial dimensions
##

conv(a::AbstractFill, b::AbstractFill) = conv(collect(a), collect(b))
conv(a::Vector, b::AbstractFill) = conv(a, collect(b))
conv(a::AbstractFill, b::Vector) = conv(collect(a), b)


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
    length(x) ≠ 1 && throw(ArgumentError("conv(::$(typeof(r)), ::$(typeof(x))) not implemented"))
    first(x)*r
end
function conv(x::AbstractVector, r::InfRanges)
    length(x) ≠ 1 && throw(ArgumentError("conv(::$(typeof(r)), ::$(typeof(x))) not implemented"))
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

