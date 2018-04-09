# This file is mmodified from Julia. License is MIT: https://julialang.org/license

colon(start::T, stop::Infinity) where {T<:Integer} = InfUnitRange{T}(start)
colon(start::T, step::T, stop::AbstractInfinity) where {T<:Real} = InfStepRange(start, step, stop)
colon(start::T, step::Real, stop::AbstractInfinity) where {T<:Real} = colon(promote(start, step)..., stop)

# AbstractFloat specializations
colon(a::T, b::AbstractInfinity) where {T<:Real} = colon(a, T(1), b)

colon(start::T, step::T, stop::Infinity) where {T<:Real} = InfStepRange(start,step,stop)

# this is needed for showarray
colon(::Infinity, ::Infinity) = 1:0

## 1-dimensional ranges ##

abstract type AbstractInfRange{T} <: InfVector{T} end


convert(::Type{T}, r::AbstractInfRange) where {T<:AbstractInfRange} =
    r isa T ? r : T(r)

## ordinal ranges

abstract type InfOrdinalRange{T,S,IS} <: AbstractInfRange{T} end
abstract type AbstractInfUnitRange{T} <: InfOrdinalRange{T,Int,Infinity} end

checkindex(::Type{Bool}, inds::AbstractInfUnitRange, i::Real) = (first(inds) <= i)

const InfIndices = Tuple{Vararg{Union{AbstractUnitRange,AbstractInfUnitRange},N}} where N
inds2string(inds::InfIndices) = join(map(string,inds), '×')


struct InfStepRange{T,S} <: InfOrdinalRange{T,S,OrientedInfinity{Bool}}
    start::T
    step::S
    stop::OrientedInfinity{Bool}
    function InfStepRange{T,S}(start::T, step::S,stop::OrientedInfinity{Bool}) where {T,S}
        sign(step) == sign(stop) || throw(ArgumentError("InfStepRange must have infinite length"))
        new{T,S}(start, step, stop)
    end
end

InfStepRange(start::T, step::S, stop::OrientedInfinity{Bool}) where {T,S} =InfStepRange{T,S}(start,step,stop)
InfStepRange{T,S}(start, step, stop) where {T,S} = InfStepRange{T,S}(convert(T,start),convert(S,step),OrientedInfinity{Bool}(stop))
InfStepRange(start, step, stop) = InfStepRange(start,step,OrientedInfinity{Bool}(stop))


struct InfUnitRange{T<:Real} <: AbstractInfUnitRange{T}
    start::T
end

"""
    OneToInf(n)

Define an `AbstractInfUnitRange` that behaves like `1:∞`, with the added
distinction that the limits are guaranteed (by the type system) to
be 1 and ∞.
"""
struct OneToInf{T<:Integer} <: AbstractInfUnitRange{T} end

OneToInf() = OneToInf{Int}()
OneTo(::Infinity) = OneToInf()

## interface implementations

size(r::AbstractInfRange) = (∞,)

isempty(r::AbstractInfRange) = false

step(r::InfStepRange) = r.step
step(r::AbstractInfUnitRange) = 1

length(r::AbstractInfRange) = ∞
unsafe_length(r::AbstractInfRange) = ∞

first(r::InfOrdinalRange{T}) where {T} = convert(T, r.start)
first(r::OneToInf{T}) where {T} = oneunit(T)

last(r::AbstractInfUnitRange) = ∞
last(r::InfStepRange) = r.stop

minimum(r::AbstractInfUnitRange) = first(r)
maximum(r::AbstractInfUnitRange) = last(r)


# Ranges are immutable
copy(r::AbstractInfRange) = r


## iteration
done(r::AbstractInfRange, i) = false

start(r::InfStepRange) = oftype(r.start + r.step, r.start)
next(r::InfStepRange{T}, i) where {T} = (convert(T,i), i+r.step)

start(r::InfUnitRange{T}) where {T} = oftype(r.start + oneunit(T), r.start)
next(r::AbstractInfUnitRange{T}, i) where {T} = (convert(T, i), i + oneunit(T))

start(r::OneToInf{T}) where {T} = oneunit(T)


## indexing

Base.offsetin(i, r::AbstractInfUnitRange) = i-first(r)
Base.nextL(L, r::AbstractInfUnitRange) = L*unsafe_length(r)
_sub2ind(inds::InfIndices, I::Integer...) = (@_inline_meta; _sub2ind_recurse(inds, 1, 1, I...))
_sub2ind(inds::Tuple{OneToInf}, i::Integer)    = i

function getindex(v::InfUnitRange{T}, i::Integer) where T
    @boundscheck i > 0 || Base.throw_boundserror(v, i)
    convert(T, first(v) + i - 1)
end
function getindex(v::OneToInf{T}, i::Integer) where T
    @boundscheck i > 0 || Base.throw_boundserror(v, i)
    convert(T, i)
end

function getindex(v::AbstractInfRange{T}, i::Integer) where T
    @boundscheck i > 0 || Base.throw_boundserror(v, i)
    convert(T, first(v) + (i - 1)*step(v))
end


getindex(r::AbstractInfRange, ::Colon) = copy(r)

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
    range(st, length(s))
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
    range(st, step(s), length(s))
end

function getindex(r::InfStepRange, s::AbstractInfRange{<:Integer})
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
    range(st, step(r)*step(s), length(s))
end

show(io::IO, r::AbstractInfRange) = print(io, repr(first(r)), ':', repr(step(r)), ':', repr(last(r)))
show(io::IO, r::AbstractInfUnitRange) = print(io, repr(first(r)), ':', repr(last(r)))
show(io::IO, r::OneToInf) = print(io, "OneToInf()")

==(r::T, s::T) where {T<:AbstractInfRange} = (first(r) == first(s)) & (step(r) == step(s))
==(r::InfOrdinalRange, s::InfOrdinalRange) = (first(r) == first(s)) & (step(r) == step(s))

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

intersect(s::StepRange, r::InfStepRange) = intersect(r, s)
intersect(s::AbstractRange, r::InfStepRange) = intersect(StepRange(s), r)
intersect(s::InfStepRange, r::AbstractRange) = intersect(s, StepRange(r))

function intersect(r1::Union{AbstractRange,AbstractInfRange}, r2::Union{AbstractRange,AbstractInfRange},
                   r3::Union{AbstractRange,AbstractInfRange}, r::Union{AbstractRange,AbstractInfRange}...)
    i = intersect(intersect(r1, r2), r3)
    for t in r
        i = intersect(i, t)
    end
    i
end

## linear operations on ranges ##
for op in (:*, :\)
    @eval $op(x::Number, r::AbstractInfRange) = $op(x,first(r)):$op(x,step(r)):$op(x,last(r))
end
for op in (:*, :/)
    @eval $op(r::AbstractInfRange, x::Number) = $op(first(r),x):$op(step(r),x):$op(last(r),x)
end

## scalar-range broadcast operations ##

broadcast(::typeof(-), r::InfOrdinalRange) = -first(r):-step(r):-last(r)

broadcast(::typeof(+), x::Real, r::AbstractInfUnitRange) = x+first(r):∞
broadcast(::typeof(+), x::Number, r::AbstractInfRange) = (x+first(r)):step(r):last(r)
broadcast(::typeof(+), r::AbstractInfRange, x::Number) = broadcast(+, x, r)  # assumes addition is commutative

broadcast(::typeof(-), x::Number, r::AbstractInfRange) = (x-first(r)):-step(r):(-last(r))
broadcast(::typeof(-), r::AbstractInfRange, x::Number) = broadcast(+, -x, r)  # assumes addition is commutative

# separate in case of noncommutative multiplication
for op in (:*, :\)
    @eval broadcast(::typeof($op), x::Number, r::AbstractInfRange) = $op(x,r)
end
for op in (:*, :/)
    @eval broadcast(::typeof($op), r::AbstractInfRange, x::Number) = $op(r,x)
end

promote_rule(a::Type{InfUnitRange{T1}}, b::Type{InfUnitRange{T2}}) where {T1,T2} =
    InfUnitRange{promote_type(T1,T2)}
InfUnitRange{T}(r::InfUnitRange{T}) where {T<:Real} = r
InfUnitRange{T}(r::UnitRange) where {T<:Real} = InfUnitRange{T}(r.start)

promote_rule(a::Type{OneToInf{T1}}, b::Type{OneToInf{T2}}) where {T1,T2} =
    OneToInf{promote_type(T1,T2)}
OneToInf{T}(r::OneToInf{T}) where {T<:Integer} = r
OneToInf{T}(r::OneToInf) where {T<:Integer} = OneToInf{T}()

promote_rule(a::Type{InfUnitRange{T1}}, ::Type{UR}) where {T1,UR<:AbstractInfUnitRange} =
    promote_rule(a, InfUnitRange{eltype(UR)})
InfUnitRange{T}(r::AbstractInfUnitRange) where {T<:Real} = InfUnitRange{T}(first(r))
InfUnitRange(r::AbstractInfUnitRange) = InfUnitRange(first(r))

AbstractInfUnitRange{T}(r::AbstractInfUnitRange{T}) where {T} = r
AbstractInfUnitRange{T}(r::InfUnitRange) where {T} = InfUnitRange{T}(r)
AbstractInfUnitRange{T}(r::OneToInf) where {T} = OneToInf{T}(r)

promote_rule(::Type{InfStepRange{T1a,T1b}}, ::Type{InfStepRange{T2a,T2b}}) where {T1a,T1b,T2a,T2b} =
    InfStepRange{promote_type(T1a,T2a), promote_type(T1b,T2b)}
InfStepRange{T1,T2}(r::InfStepRange{T1,T2}) where {T1,T2} = r

promote_rule(a::Type{InfStepRange{T1a,T1b}}, ::Type{UR}) where {T1a,T1b,UR<:AbstractInfUnitRange} =
    promote_rule(a, InfStepRange{eltype(UR), eltype(UR)})
InfStepRange{T1,T2}(r::AbstractInfRange) where {T1,T2} =
    InfStepRange{T1,T2}(convert(T1, first(r)), convert(T2, step(r)), last(r))
InfStepRange(r::AbstractInfUnitRange{T}) where {T} =
    InfStepRange{T,T}(first(r), step(r), last(r))
(::Type{InfStepRange{T1,T2} where T1})(r::AbstractInfRange) where {T2} =
    InfStepRange{eltype(r),T2}(r)





## sorting ##

issorted(r::AbstractInfUnitRange) = true
issorted(r::AbstractInfRange) = true

sort(r::AbstractInfUnitRange) = r
sort!(r::AbstractInfUnitRange) = r

sort(r::AbstractInfRange) = issorted(r) ? r : throw(ArgumentError("Cannot sort infinite range"))

sortperm(r::AbstractInfUnitRange) = 1:∞
sortperm(r::AbstractInfRange) = 1:1:∞

sum(r::AbstractInfRange{<:Real}) = last(r)

mean(r::AbstractInfRange{<:Real}) = last(r)

median(r::AbstractInfRange{<:Real}) = last(r)

function _in_range(x, r::AbstractInfRange)
    if step(r) == 0
        return !isempty(r) && first(r) == x
    else
        n = round(Integer, (x - first(r)) / step(r)) + 1
        return n >= 1 && r[n] == x
    end
end
in(x::AbstractInfinity, r::AbstractInfRange) = false # never reach it...
in(x::Real, r::AbstractInfRange{<:Real}) = _in_range(x, r)
# This method needs to be defined separately since -(::T, ::T) can be implemented
# even if -(::T, ::Real) is not
in(x::T, r::AbstractInfRange{T}) where {T} = _in_range(x, r)

in(x::Integer, r::AbstractInfUnitRange{<:Integer}) = (first(r) <= x)

in(x::Real, r::AbstractInfRange{T}) where {T<:Integer} =
    isinteger(x) && !isempty(r) && ifelse(step(r) > 0, x ≥ first(r), x ≤ first(r)) &&
            (mod(convert(T,x),step(r))-mod(first(r),step(r)) == 0)
# Addition/subtraction of ranges

-(r1::AbstractInfUnitRange, r2::AbstractInfUnitRange) = Fill(first(r1)-first(r2), ∞)

for op in (:+, :-)
    @eval begin
        function $op(r1::InfOrdinalRange, r2::InfOrdinalRange)
            new_step = $op(step(r1), step(r2))
            $op(first(r1), first(r2)):new_step:sign(new_step)*∞
        end
        broadcast(::typeof($op), r1::AbstractInfRange, r2::AbstractInfRange) = $op(r1, r2)
    end
end
