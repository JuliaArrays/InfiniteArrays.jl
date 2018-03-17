# This file is mmodified from Julia. License is MIT: https://julialang.org/license

colon(start::T, stop::Infinity) where {T<:Integer} = InfiniteUnitRange{T}(start)
colon(start::T, step::T, stop::AbstractInfinity) where {T<:Real} = InfiniteStepRange(start, step, stop)
colon(start::T, step::Real, stop::AbstractInfinity) where {T<:Real} = colon(promote(start, step)..., stop)

# AbstractFloat specializations
colon(a::T, b::AbstractInfinity) where {T<:Real} = colon(a, T(1), b)

colon(start::T, step::T, stop::Infinity) where {T<:Real} = InfiniteStepRange(start,step,stop)

## 1-dimensional ranges ##

abstract type AbstractInfiniteRange{T} <: InfiniteVector{T} end


convert(::Type{T}, r::AbstractInfiniteRange) where {T<:AbstractInfiniteRange} =
    r isa T ? r : T(r)

## ordinal ranges

abstract type InfiniteOrdinalRange{T,S,IS} <: AbstractInfiniteRange{T} end
abstract type AbstractInfiniteUnitRange{T} <: InfiniteOrdinalRange{T,Int,Infinity} end

struct InfiniteStepRange{T,S} <: InfiniteOrdinalRange{T,S,OrientedInfinity{Bool}}
    start::T
    step::S
    stop::OrientedInfinity{Bool}
    function InfiniteStepRange{T,S}(start::T, step::S,stop::OrientedInfinity{Bool}) where {T,S}
        sign(step) == sign(stop) || throw(ArgumentError("InfiniteStepRange must have infinite length"))
        new{T,S}(start, step, stop)
    end
end

InfiniteStepRange(start::T, step::S, stop::OrientedInfinity{Bool}) where {T,S} =InfiniteStepRange{T,S}(start,step,stop)
InfiniteStepRange{T,S}(start, step, stop) where {T,S} = InfiniteStepRange{T,S}(convert(T,start),convert(S,step),OrientedInfinity{Bool}(stop))
InfiniteStepRange(start, step, stop) = InfiniteStepRange(start,step,OrientedInfinity{Bool}(stop))


struct InfiniteUnitRange{T<:Real} <: AbstractInfiniteUnitRange{T}
    start::T
end

"""
    OneToInfinity(n)

Define an `AbstractInfiniteUnitRange` that behaves like `1:∞`, with the added
distinction that the limits are guaranteed (by the type system) to
be 1 and ∞.
"""
struct OneToInfinity{T<:Integer} <: AbstractInfiniteUnitRange{T} end

OneToInfinity() = OneToInfinity{Int}()

## interface implementations

size(r::AbstractInfiniteRange) = (∞,)

isempty(r::AbstractInfiniteRange) = false

step(r::InfiniteStepRange) = r.step
step(r::AbstractInfiniteUnitRange) = 1

length(r::AbstractInfiniteRange) = ∞


first(r::InfiniteOrdinalRange{T}) where {T} = convert(T, r.start)
first(r::OneToInfinity{T}) where {T} = oneunit(T)

last(r::AbstractInfiniteUnitRange) = ∞
last(r::InfiniteStepRange) = r.stop

minimum(r::AbstractInfiniteUnitRange) = first(r)
maximum(r::AbstractInfiniteUnitRange) = last(r)


# Ranges are immutable
copy(r::AbstractInfiniteRange) = r


## iteration
done(r::AbstractInfiniteRange, i) = false

start(r::InfiniteStepRange) = oftype(r.start + r.step, r.start)
next(r::InfiniteStepRange{T}, i) where {T} = (convert(T,i), i+r.step)

start(r::InfiniteUnitRange{T}) where {T} = oftype(r.start + oneunit(T), r.start)
next(r::AbstractInfiniteUnitRange{T}, i) where {T} = (convert(T, i), i + oneunit(T))

start(r::OneToInfinity{T}) where {T} = oneunit(T)


## indexing

function getindex(v::InfiniteUnitRange{T}, i::Integer) where T
    @boundscheck i > 0 || Base.throw_boundserror(v, i)
    convert(T, first(v) + i - 1)
end
function getindex(v::OneToInfinity{T}, i::Integer) where T
    @boundscheck i > 0 || Base.throw_boundserror(v, i)
    convert(T, i)
end

function getindex(v::AbstractInfiniteRange{T}, i::Integer) where T
    @boundscheck i > 0 || Base.throw_boundserror(v, i)
    convert(T, first(v) + (i - 1)*step(v))
end


getindex(r::AbstractInfiniteRange, ::Colon) = copy(r)

function getindex(r::AbstractInfiniteUnitRange, s::AbstractInfiniteUnitRange{<:Integer})
    f = first(r)
    @boundscheck first(s) ≥ 1 || Base.throw_boundserror(r, first(s))
    st = oftype(f, f + first(s)-1)
    st:∞
end

function getindex(r::AbstractInfiniteUnitRange, s::AbstractUnitRange{<:Integer})
    f = first(r)
    @boundscheck first(s) ≥ 1 || Base.throw_boundserror(r, first(s))
    st = oftype(f, f + first(s)-1)
    range(st, length(s))
end

getindex(r::OneToInfinity{T}, s::OneTo) where T = OneTo(T(s.stop))

function getindex(r::AbstractInfiniteUnitRange, s::InfiniteStepRange{<:Integer})
    @_inline_meta
    @boundscheck (step(s) > 0 && first(s) ≥ 1) || throw(BoundsError(r, minimum(s)))
    st = oftype(first(r), first(r) + s.start-1)
    new_step = step(s)
    st:new_step:sign(new_step)*∞
end

function getindex(r::AbstractInfiniteUnitRange, s::StepRange{<:Integer})
    @_inline_meta
    @boundscheck minimum(s) ≥ 1 || throw(BoundsError(r, minimum(s)))
    st = oftype(first(r), first(r) + s.start-1)
    range(st, step(s), length(s))
end

function getindex(r::InfiniteStepRange, s::AbstractInfiniteRange{<:Integer})
    @_inline_meta
    @boundscheck (step(s) > 0 && first(s) ≥ 1) || throw(BoundsError(r, minimum(s)))
    st = oftype(r.start, r.start + (first(s)-1)*step(r))
    new_step = step(r)*step(s)
    st:new_step:sign(new_step)*∞
end

function getindex(r::InfiniteStepRange, s::AbstractRange{<:Integer})
    @_inline_meta
    @boundscheck isempty(s) || minimum(s) ≥ 1 || throw(BoundsError(r, minimum(s)))
    st = oftype(r.start, r.start + (first(s)-1)*step(r))
    range(st, step(r)*step(s), length(s))
end

show(io::IO, r::AbstractInfiniteRange) = print(io, repr(first(r)), ':', repr(step(r)), ':', repr(last(r)))
show(io::IO, r::AbstractInfiniteUnitRange) = print(io, repr(first(r)), ':', repr(last(r)))
show(io::IO, r::OneToInfinity) = print(io, "OneToInfinity()")

==(r::T, s::T) where {T<:AbstractInfiniteRange} = (first(r) == first(s)) & (step(r) == step(s))
==(r::InfiniteOrdinalRange, s::InfiniteOrdinalRange) = (first(r) == first(s)) & (step(r) == step(s))

intersect(r::OneToInfinity{T}, s::OneToInfinity{V}) where {T,V} = OneToInfinity{promote_type(T,V)}()
intersect(r::OneToInfinity{T}, s::OneTo{T}) where T = s
intersect(r::OneTo{T}, s::OneToInfinity{T}) where T = r

intersect(r::AbstractInfiniteUnitRange{<:Integer}, s::AbstractInfiniteUnitRange{<:Integer}) =
    InfiniteUnitRange(max(first(r),first(s)))
intersect(r::AbstractInfiniteUnitRange{<:Integer}, s::AbstractUnitRange{<:Integer}) =
    max(first(r),first(s)):last(s)
intersect(r::AbstractUnitRange{<:Integer}, s::AbstractInfiniteUnitRange{<:Integer}) =
    intersect(s, r)


intersect(i::Integer, r::AbstractInfiniteUnitRange{<:Integer}) =
    i < first(r) ? (first(r):i) : (i:i)

intersect(r::AbstractInfiniteUnitRange{<:Integer}, i::Integer) = intersect(i, r)

function intersect(r::AbstractInfiniteUnitRange{<:Integer}, s::StepRange{<:Integer})
    if isempty(s)
        range(first(r), 0)
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

intersect(r::StepRange{<:Integer}, s::AbstractInfiniteUnitRange{<:Integer}) =
    intersect(s, r)


function intersect(r::InfiniteStepRange, s::InfiniteStepRange)
    sign(step(r)) == sign(step(s)) || throw(ArgumentError("Can only take intersection of infinite ranges"))
    start1 = first(r)
    step1 = step(r)
    start2 = first(s)
    step2 = step(s)
    a = lcm(step1, step2)

    g, x, y = gcdx(step1, step2)

    if rem(start1 - start2, g) != 0
        # Unaligned, no overlap possible.
        throw(ArgumentError("Cannot take intersection of InfiniteStepRange that has no elements"))
    end

    z = div(start1 - start2, g)
    b = start1 - x * z * step1
    # Possible points of the intersection of r and s are
    # ..., b-2a, b-a, b, b+a, b+2a, ...
    # Determine where in the sequence to start and stop.
    m = max(start1 + mod(b - start1, a), start2 + mod(b - start2, a))
    m:a:∞
end

intersect(r::InfiniteStepRange, s::AbstractInfiniteUnitRange) = intersect(r, InfiniteStepRange(s))
intersect(r::AbstractInfiniteUnitRange, s::InfiniteStepRange) = intersect(InfiniteStepRange(r), s)

function intersect(r::InfiniteStepRange, s::StepRange)
    if isempty(s)
        return range(first(r), step(r), 0)
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
        return range(start1, a, 0)
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

intersect(s::StepRange, r::InfiniteStepRange) = intersect(r, s)
intersect(s::AbstractRange, r::InfiniteStepRange) = intersect(StepRange(s), r)
intersect(s::InfiniteStepRange, r::AbstractRange) = intersect(s, StepRange(r))

function intersect(r1::Union{AbstractRange,AbstractInfiniteRange}, r2::Union{AbstractRange,AbstractInfiniteRange},
                   r3::Union{AbstractRange,AbstractInfiniteRange}, r::Union{AbstractRange,AbstractInfiniteRange}...)
    i = intersect(intersect(r1, r2), r3)
    for t in r
        i = intersect(i, t)
    end
    i
end

## linear operations on ranges ##
for op in (:*, :\)
    @eval $op(x::Number, r::AbstractInfiniteRange) = $op(x,first(r)):$op(x,step(r)):$op(x,last(r))
end
for op in (:*, :/)
    @eval $op(r::AbstractInfiniteRange, x::Number) = $op(first(r),x):$op(step(r),x):$op(last(r),x)
end

## scalar-range broadcast operations ##

broadcast(::typeof(-), r::InfiniteOrdinalRange) = -first(r):-step(r):-last(r)

broadcast(::typeof(+), x::Real, r::AbstractInfiniteUnitRange) = x+first(r):∞
broadcast(::typeof(+), x::Number, r::AbstractInfiniteRange) = (x+first(r)):step(r):last(r)
broadcast(::typeof(+), r::AbstractInfiniteRange, x::Number) = broadcast(+, x, r)  # assumes addition is commutative

broadcast(::typeof(-), x::Number, r::AbstractInfiniteRange) = (x-first(r)):-step(r):(-last(r))
broadcast(::typeof(-), r::AbstractInfiniteRange, x::Number) = broadcast(+, -x, r)  # assumes addition is commutative

# separate in case of noncommutative multiplication
for op in (:*, :\)
    @eval broadcast(::typeof($op), x::Number, r::AbstractInfiniteRange) = $op(x,r)
end
for op in (:*, :/)
    @eval broadcast(::typeof($op), r::AbstractInfiniteRange, x::Number) = $op(r,x)
end

promote_rule(a::Type{InfiniteUnitRange{T1}}, b::Type{InfiniteUnitRange{T2}}) where {T1,T2} =
    InfiniteUnitRange{promote_type(T1,T2)}
InfiniteUnitRange{T}(r::InfiniteUnitRange{T}) where {T<:Real} = r
InfiniteUnitRange{T}(r::UnitRange) where {T<:Real} = InfiniteUnitRange{T}(r.start)

promote_rule(a::Type{OneToInfinity{T1}}, b::Type{OneToInfinity{T2}}) where {T1,T2} =
    OneToInfinity{promote_type(T1,T2)}
OneToInfinity{T}(r::OneToInfinity{T}) where {T<:Integer} = r
OneToInfinity{T}(r::OneToInfinity) where {T<:Integer} = OneToInfinity{T}()

promote_rule(a::Type{InfiniteUnitRange{T1}}, ::Type{UR}) where {T1,UR<:AbstractInfiniteUnitRange} =
    promote_rule(a, InfiniteUnitRange{eltype(UR)})
InfiniteUnitRange{T}(r::AbstractInfiniteUnitRange) where {T<:Real} = InfiniteUnitRange{T}(first(r))
InfiniteUnitRange(r::AbstractInfiniteUnitRange) = InfiniteUnitRange(first(r))

AbstractInfiniteUnitRange{T}(r::AbstractInfiniteUnitRange{T}) where {T} = r
AbstractInfiniteUnitRange{T}(r::InfiniteUnitRange) where {T} = InfiniteUnitRange{T}(r)
AbstractInfiniteUnitRange{T}(r::OneToInfinity) where {T} = OneToInfinity{T}(r)

promote_rule(::Type{InfiniteStepRange{T1a,T1b}}, ::Type{InfiniteStepRange{T2a,T2b}}) where {T1a,T1b,T2a,T2b} =
    InfiniteStepRange{promote_type(T1a,T2a), promote_type(T1b,T2b)}
InfiniteStepRange{T1,T2}(r::InfiniteStepRange{T1,T2}) where {T1,T2} = r

promote_rule(a::Type{InfiniteStepRange{T1a,T1b}}, ::Type{UR}) where {T1a,T1b,UR<:AbstractInfiniteUnitRange} =
    promote_rule(a, InfiniteStepRange{eltype(UR), eltype(UR)})
InfiniteStepRange{T1,T2}(r::AbstractInfiniteRange) where {T1,T2} =
    InfiniteStepRange{T1,T2}(convert(T1, first(r)), convert(T2, step(r)), last(r))
InfiniteStepRange(r::AbstractInfiniteUnitRange{T}) where {T} =
    InfiniteStepRange{T,T}(first(r), step(r), last(r))
(::Type{InfiniteStepRange{T1,T2} where T1})(r::AbstractInfiniteRange) where {T2} =
    InfiniteStepRange{eltype(r),T2}(r)


## sorting ##

issorted(r::AbstractInfiniteUnitRange) = true
issorted(r::AbstractInfiniteRange) = true

sort(r::AbstractInfiniteUnitRange) = r
sort!(r::AbstractInfiniteUnitRange) = r

sort(r::AbstractInfiniteRange) = issorted(r) ? r : throw(ArgumentError("Cannot sort infinite range"))

sortperm(r::AbstractInfiniteUnitRange) = 1:∞
sortperm(r::AbstractInfiniteRange) = 1:1:∞

sum(r::AbstractInfiniteRange{<:Real}) = last(r)

mean(r::AbstractInfiniteRange{<:Real}) = last(r)

median(r::AbstractInfiniteRange{<:Real}) = last(r)

function _in_range(x, r::AbstractInfiniteRange)
    if step(r) == 0
        return !isempty(r) && first(r) == x
    else
        n = round(Integer, (x - first(r)) / step(r)) + 1
        return n >= 1 && r[n] == x
    end
end
in(x::AbstractInfinity, r::AbstractInfiniteRange) = false # never reach it...
in(x::Real, r::AbstractInfiniteRange{<:Real}) = _in_range(x, r)
# This method needs to be defined separately since -(::T, ::T) can be implemented
# even if -(::T, ::Real) is not
in(x::T, r::AbstractInfiniteRange{T}) where {T} = _in_range(x, r)

in(x::Integer, r::AbstractInfiniteUnitRange{<:Integer}) = (first(r) <= x)

in(x::Real, r::AbstractInfiniteRange{T}) where {T<:Integer} =
    isinteger(x) && !isempty(r) && ifelse(step(r) > 0, x ≥ first(r), x ≤ first(r)) &&
            (mod(convert(T,x),step(r))-mod(first(r),step(r)) == 0)
# Addition/subtraction of ranges

-(r1::AbstractInfiniteUnitRange, r2::AbstractInfiniteUnitRange) = repeated(first(r1)-first(r2))

for op in (:+, :-)
    @eval begin
        function $op(r1::InfiniteOrdinalRange, r2::InfiniteOrdinalRange)
            new_step = $op(step(r1), step(r2))
            $op(first(r1), first(r2)):new_step:sign(new_step)*∞
        end
        broadcast(::typeof($op), r1::AbstractInfiniteRange, r2::AbstractInfiniteRange) = $op(r1, r2)
    end
end
