# This file is mmodified from Julia. License is MIT: https://julialang.org/license

const PosInfinity = Union{Infinity, InfiniteCardinal{0}}

(:)(start::T, stop::PosInfinity) where {T<:Integer} = InfUnitRange{T}(start)
(:)(start::PosInfinity, stop::Integer) = start+1:start
function (:)(start::T, step::T, stop::RealInfinity) where {T<:Real}
    signbit(step) == signbit(stop) || throw(ArgumentError("InfStepRange must have infinite length"))
    InfStepRange(start, step)
end
(:)(start::T, step::Real, stop::RealInfinity) where {T<:Real} = (:)(promote(start, step)..., stop)
(:)(start::Real, step, stop::PosInfinity)= (:)(start, step, RealInfinity(stop))
(:)(::PosInfinity, _, ::Real) = throw(ArgumentError("Cannot create range starting at infinity"))

# AbstractFloat specializations
(:)(a::T, b::Union{PosInfinity,RealInfinity}) where {T<:Real} = (:)(a, T(1), b)

function (:)(start::T, step::T, stop::PosInfinity) where {T<:Real}
    sign(step) == sign(stop) || throw(ArgumentError("InfStepRange must have infinite length"))
    InfStepRange(start,step)
end

# this is needed for showarray
(:)(::PosInfinity, ::PosInfinity) = 1:0


# Range of a given length: range(a, [step=s,] length=l), no stop
_range(a::Real,          ::Nothing,         ::Nothing, len::InfiniteCardinal{0}) = InfUnitRange{typeof(a)}(a)
_range(a::AbstractFloat, ::Nothing,         ::Nothing, len::InfiniteCardinal{0}) = _range(a, oftype(a, 1),   nothing, len)
_range(a::T, st::T, ::Nothing, ::InfiniteCardinal{0}) where T<:IEEEFloat = InfStepRange{T,T}(a, st)
_range(a::T, st::T, ::Nothing, ::InfiniteCardinal{0}) where T<:AbstractFloat = InfStepRange{T,T}(a, st)
if VERSION < v"1.8-"
    _rangestyle(::Ordered, ::ArithmeticWraps, a::T, step::S, len::InfiniteCardinal{0}) where {T,S} = InfStepRange{T,S}(a, step)
    _rangestyle(::Ordered, ::ArithmeticUnknown, a::T, step::S, len::InfiniteCardinal{0}) where {T,S} = InfStepRange{T,S}(a, step)
    _rangestyle(::Any, ::Any, a::T, step::S, len::InfiniteCardinal{0}) where {T,S} = InfStepRange{T,S}(a, step)
else
    range_start_step_length(a, st, ::InfiniteCardinal{0}) = InfStepRange(a, st)
    range_start_step_length(a::Real, st::IEEEFloat, len::InfiniteCardinal{0}) = range_start_step_length(promote(a, st)..., len)
    range_start_step_length(a::IEEEFloat, st::Real, len::InfiniteCardinal{0}) = range_start_step_length(promote(a, st)..., len)
    range_start_step_length(a::IEEEFloat, st::IEEEFloat, len::InfiniteCardinal{0}) = range_start_step_length(promote(a, st)..., len)
    range_start_step_length(a::T, st::T, ::InfiniteCardinal{0}) where T<:IEEEFloat = InfStepRange{T,T}(a, st)
end


# Construct range for rational start=start_n/den, step=step_n/den
floatrange(::Type{T}, start_n::Integer, step_n::Integer, ::PosInfinity, den::Integer) where T =
    InfStepRange(T(start_n)/den,T(step_n)/den)

floatrange(a::AbstractFloat, st::AbstractFloat, ::PosInfinity, divisor::AbstractFloat) =
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

FillArrays.steprangelen(start, step, ::PosInfinity) = InfStepRange(start, step)

abstract type AbstractInfUnitRange{T<:Real} <: AbstractUnitRange{T} end

done(r::AbstractInfUnitRange{T}, i) where {T} = false
unitrange_last(start, stop::PosInfinity) = ∞

struct InfUnitRange{T<:Real} <: AbstractInfUnitRange{T}
    start::T
end


InfUnitRange(a::InfUnitRange) = a
InfUnitRange{T}(a::AbstractInfUnitRange) where T<:Real = InfUnitRange{T}(first(a))
InfUnitRange(a::AbstractInfUnitRange{T}) where T<:Real = InfUnitRange{T}(first(a))
unitrange(a::AbstractInfUnitRange) = InfUnitRange(a)

AbstractArray{T}(a::InfUnitRange) where T<:Real = InfUnitRange{T}(a.start)
AbstractVector{T}(a::InfUnitRange) where T<:Real = InfUnitRange{T}(a.start)
AbstractArray{T}(a::InfStepRange) where T<:Real = InfStepRange(convert(T,a.start), convert(T,a.step))
AbstractVector{T}(a::InfStepRange) where T<:Real = InfStepRange(convert(T,a.start), convert(T,a.step))

const InfRanges{T} = Union{InfStepRange{T},AbstractInfUnitRange{T}}
const InfAxes = Union{InfRanges{<:Integer},Slice{<:AbstractInfUnitRange{<:Integer}},IdentityUnitRange{<:AbstractInfUnitRange{<:Integer}}}

Base.IteratorSize(::Type{<:InfRanges}) = Base.IsInfinite()

AbstractArray{T}(ac::Adjoint{<:Any,<:InfRanges}) where T<:Real = AbstractArray{T}(parent(ac))'
AbstractMatrix{T}(ac::Adjoint{<:Any,<:InfRanges}) where T<:Real = AbstractVector{T}(parent(ac))'
AbstractArray{T}(ac::Transpose{<:Any,<:InfRanges}) where T<:Real = transpose(AbstractArray{T}(parent(ac)))
AbstractMatrix{T}(ac::Transpose{<:Any,<:InfRanges}) where T<:Real = transpose(AbstractVector{T}(parent(ac)))

copy(ac::AdjOrTrans{<:Any,<:InfRanges}) = ac
reverse(a::InfRanges) = throw(ArgumentError("Cannot reverse infinite range"))

if VERSION ≥ v"1.7-"
    Base.range_start_step_length(a::T, st::T, len::PosInfinity) where T<:Union{Float16,Float32,Float64} = InfStepRange(a, st)
    for op in (:+, :-)
        @eval $op(a::InfRanges, b::InfRanges) = InfStepRange($op(first(a),first(b)), $op(step(a),step(b)))
    end
    broadcasted(::DefaultArrayStyle{1}, ::typeof(*), x::Number, r::InfRanges) = InfStepRange(x*first(r), x*step(r))
    broadcasted(::DefaultArrayStyle{1}, ::typeof(*), r::InfRanges, x::Number) = InfStepRange(first(r)*x, step(r)*x)
    broadcasted(::DefaultArrayStyle{1}, ::typeof(*), x::AbstractFloat, r::InfRanges) = InfStepRange(x*first(r), x*step(r))
    broadcasted(::DefaultArrayStyle{1}, ::typeof(*), r::InfRanges, x::AbstractFloat) = InfStepRange(first(r)*x, step(r)*x)
end


"""
    OneToInf(n)

Define an `AbstractInfUnitRange` that behaves like `1:∞`, with the added
distinction that the limits are guaranteed (by the type system) to
be 1 and ∞.
"""
struct OneToInf{T<:Integer} <: AbstractInfUnitRange{T} end

OneToInf() = OneToInf{Int}()
oneto(::PosInfinity) = OneToInf()
function oneto(x::ComplexInfinity)
    iszero(angle(x)) && return oneto(∞)
    throw(ArgumentError("Cannot create infinite range with negative length"))
 end
 function oneto(x::RealInfinity)
    signbit(x) || return oneto(∞)
    throw(ArgumentError("Cannot create infinite range with negative length"))
 end

AbstractArray{T}(a::OneToInf) where T<:Integer = OneToInf{T}()
AbstractVector{T}(a::OneToInf) where T<:Integer = OneToInf{T}()
AbstractArray{T}(a::OneToInf) where T<:Real = InfUnitRange{T}(a)
AbstractVector{T}(a::OneToInf) where T<:Real = InfUnitRange{T}(a)


(==)(::OneToInf, ::OneToInf) = true

## interface implementations

size(r::InfRanges) = (ℵ₀,)
axes(r::InfRanges{<:Integer}) = (OneToInf{promote_type(Int,eltype(r))}(),)
axes(r::InfRanges) = (OneToInf(),)

isempty(r::InfRanges) = false

step(r::InfStepRange) = r.step

length(r::InfRanges) = ℵ₀
if VERSION >= v"1.7"
    Base.checked_length(r::InfRanges) = ℵ₀
end
unsafe_length(r::InfRanges) = ℵ₀

first(r::OneToInf{T}) where {T} = oneunit(T)

last(r::AbstractInfUnitRange) = ℵ₀
last(r::InfStepRange) = RealInfinity(signbit(step(r)))

minimum(r::InfUnitRange) = first(r)
maximum(r::InfUnitRange) = ℵ₀


## iteration
start(r::InfStepRange) = oftype(r.start + r.step, r.start)
next(r::InfStepRange{T}, i) where {T} = (convert(T,i), i+r.step)

start(r::InfUnitRange{T}) where {T} = oftype(r.start + oneunit(T), r.start)
start(r::OneToInf{T}) where {T} = oneunit(T)

done(r::InfStepRange{T}, i) where {T} = false

## iteration with zip + finite iterator
#  allows axes(zip(...)) and size(zip(...))
Base.Iterators._promote_tuple_shape((a,)::Tuple{OneToInf}, (b,)::Tuple{OneTo}) =
    (intersect(a, b),)
Base.Iterators._promote_tuple_shape((a,)::Tuple{OneTo}, (b,)::Tuple{OneToInf}) =
    (intersect(a, b),)

## indexing

unsafe_indices(S::Slice{<:OneToInf}) = (S.indices,)
axes(S::Slice{<:OneToInf}) = (S.indices,)
_sub2ind(inds::Tuple{OneToInf}, i::Integer)    = i

to_shape(::OneToInf) = ℵ₀

# used for linear indexing
function _ind2sub_recurse(inds::Tuple{OneToInf{Int},Vararg{Any}}, ind::Integer)
    @_inline_meta
    (ind+1, _ind2sub_recurse(tail(inds), 0)...)
end

function _ind2sub_recurse(indslast::Tuple{OneToInf{Int}}, ind::Integer)
	@_inline_meta
	(ind+1,)
end

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
    i == 1 && return convert(T, first(v))
    convert(T, first(v) + (i - 1)*step(v))
end

function getindex(x::AbstractUnitRange, y::PosInfinity)
    isinf(length(x)) || throw(BoundsError(x,y))
    ℵ₀
end

function getindex(x::OneToInf{T}, y::PosInfinity) where T
    isinf(length(x)) || throw(BoundsError(x,y))
    ℵ₀
end

function getindex(x::InfStepRange{T}, y::PosInfinity) where T
    isinf(length(x)) || throw(BoundsError(x,y))
    ℵ₀
end
function getindex(x::InfUnitRange{T}, y::PosInfinity) where T
    isinf(length(x)) || throw(BoundsError(x,y))
    ℵ₀
end

getindex(::AbstractInfUnitRange, ::Infinity) = ℵ₀
getindex(::OneToInf, ::Infinity) = ℵ₀
getindex(v::InfUnitRange, i::Infinity) = ℵ₀
getindex(v::InfStepRange, i::Infinity) = ℵ₀

function getindex(r::AbstractInfUnitRange, s::AbstractInfUnitRange{<:Integer})
    f = first(r)
    @boundscheck first(s) ≥ 1 || Base.throw_boundserror(r, first(s))
    st = oftype(f, f + first(s)-1)
    st:ℵ₀
end

function getindex(r::AbstractInfUnitRange, s::AbstractUnitRange{<:Integer})
    f = first(r)
    @boundscheck first(s) ≥ 1 || Base.throw_boundserror(r, first(s))
    st = oftype(f, f + first(s)-1)
    range(st; length=length(s))
end

getindex(r::AbstractInfUnitRange, s::Slice{<:AbstractInfUnitRange{<:Integer}}) = r

getindex(r::OneToInf{T}, s::OneTo) where T = OneTo(T(s.stop))

function getindex(r::AbstractInfUnitRange, s::InfStepRange{<:Integer})
    @_inline_meta
    @boundscheck (step(s) > 0 && first(s) ≥ 1) || throw(BoundsError(r, minimum(s)))
    st = oftype(first(r), first(r) + s.start-1)
    new_step = step(s)
    InfStepRange(st,new_step)
end

function getindex(r::AbstractInfUnitRange, s::StepRange{<:Integer})
    @_inline_meta
    @boundscheck minimum(s) ≥ 1 || throw(BoundsError(r, minimum(s)))
    st = oftype(first(r), first(r) + s.start-1)
    range(st; step=step(s), length=length(s))
end

function getindex(r::InfStepRange, s::InfAxes)
    @_inline_meta
    @boundscheck (step(s) > 0 && first(s) ≥ 1) || throw(BoundsError(r, minimum(s)))
    st = oftype(r.start, r.start + (first(s)-1)*step(r))
    new_step = step(r)*step(s)
    InfStepRange(st,new_step)
end

function getindex(r::InfStepRange, s::AbstractRange{<:Integer})
    @_inline_meta
    @boundscheck isempty(s) || minimum(s) ≥ 1 || throw(BoundsError(r, minimum(s)))
    st = oftype(r.start, r.start + (first(s)-1)*step(r))
    range(st; step=step(r)*step(s), length=length(s))
end

function getindex(Ac::AdjOrTrans{<:Any,<:InfRanges}, k::Integer, j)
    @boundscheck k == 1 || throw(BoundsError(Ac, k))
    parent(Ac)[j]
end

function getindex(Ac::AdjOrTrans{<:Any,<:InfRanges}, k::Integer, j::InfAxes)
    @boundscheck k == 1 || throw(BoundsError(Ac, k))
    parent(Ac)[j]
end


show(io::IO, r::InfUnitRange) = print(io, repr(first(r)), ':', repr(∞))
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
    sign(step(r)) == sign(step(s)) || throw(ArgumentError("Can only take intersection of infinite ranges"))
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

function union(ain::InfRanges, bin::InfRanges)
    a,b = promote(ain, bin)
    first(b) ∈ a && iszero(mod(step(b), step(a))) && return a
    first(a) ∈ b && iszero(mod(step(a), step(b))) && return b
    throw(ArgumentError("Cannot take union of $a and $b"))
end

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

in(x::Union{Infinity,RealInfinity}, r::InfRanges) = false # never reach it...
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

function sort!(a::InfStepRange)
    step(a) > 0 || throw(ArgumentError("Cannot sort $a"))
    a
end


### Lazy broadcasting

BroadcastStyle(::Type{<:InfRanges}) = LazyArrayStyle{1}()
BroadcastStyle(::Type{<:Adjoint{<:Any,<:InfRanges}}) = LazyArrayStyle{2}()
BroadcastStyle(::Type{<:Transpose{<:Any,<:InfRanges}}) = LazyArrayStyle{2}()
const InfSubVector = SubArray{<:Any,1,<:Any,<:Tuple{InfRanges}}
BroadcastStyle(::Type{<:InfSubVector}) = LazyArrayStyle{1}()
BroadcastStyle(::Type{<:SubArray{<:Any,1,<:Any,<:Tuple{InfSubVector}}}) = LazyArrayStyle{1}()

BroadcastStyle(::Type{<:Base.Slice{<:InfRanges}}) = LazyArrayStyle{1}()


const InfIndexRanges{T<:Integer} = Union{InfStepRange{T},AbstractInfUnitRange{T},Slice{OneToInf{T}}}

BroadcastStyle(::Type{<:SubArray{<:Any,1,<:Any,Tuple{<:InfIndexRanges}}})= LazyArrayStyle{1}()
BroadcastStyle(::Type{<:SubArray{<:Any,1,<:Any,<:Tuple{<:InfIndexRanges,<:Any}}})= LazyArrayStyle{1}()
BroadcastStyle(::Type{<:SubArray{<:Any,1,<:Any,<:Tuple{<:Any,<:InfIndexRanges}}})= LazyArrayStyle{1}()
BroadcastStyle(::Type{<:SubArray{<:Any,1,<:LazyMatrix,<:Tuple{<:InfIndexRanges,<:Any}}})= LazyArrayStyle{1}()
BroadcastStyle(::Type{<:SubArray{<:Any,1,<:LazyMatrix,<:Tuple{<:Any,<:InfIndexRanges}}})= LazyArrayStyle{1}()
BroadcastStyle(::Type{<:SubArray{<:Any,2,<:Any,<:Tuple{<:InfIndexRanges,<:InfIndexRanges}}})= LazyArrayStyle{2}()
BroadcastStyle(::Type{<:SubArray{<:Any,2,<:Any,<:Tuple{<:InfIndexRanges,<:Any}}})= LazyArrayStyle{2}()
BroadcastStyle(::Type{<:SubArray{<:Any,2,<:Any,<:Tuple{<:Any,<:InfIndexRanges}}})= LazyArrayStyle{2}()




broadcasted(::BroadcastStyle, f, r::Adjoint{<:Any,<:InfRanges}) = broadcast(f,parent(r))'
broadcasted(::BroadcastStyle, f, r::Transpose{<:Any,<:InfRanges}) = transpose(broadcast(f,parent(r)))
broadcasted(::BroadcastStyle, f, a::Number, r::Adjoint{<:Any,<:InfRanges}) = broadcast(f,a,parent(r))'
broadcasted(::BroadcastStyle, f, a::Number, r::Transpose{<:Any,<:InfRanges}) = transpose(broadcast(f,a,parent(r)))
broadcasted(::BroadcastStyle, f, r::Adjoint{<:Any,<:InfRanges}, a::Number) = broadcast(f,parent(r),a)'
broadcasted(::BroadcastStyle, f, r::Transpose{<:Any,<:InfRanges}, a::Number) = transpose(broadcast(f,parent(r),a))

broadcast(f, r::Adjoint{<:Any,<:InfRanges}) = broadcast(f,parent(r))'
broadcast(f, r::Transpose{<:Any,<:InfRanges}) = transpose(broadcast(f,parent(r)))
broadcast(f, a::Number, r::Adjoint{<:Any,<:InfRanges}) = broadcast(f,a,parent(r))'
broadcast(f, a::Number, r::Transpose{<:Any,<:InfRanges}) = transpose(broadcast(f,a,parent(r)))
broadcast(f, r::Adjoint{<:Any,<:InfRanges}, a::Number) = broadcast(f,parent(r),a)'
broadcast(f, r::Transpose{<:Any,<:InfRanges}, a::Number) = transpose(broadcast(f,parent(r),a))


# cumsum(r::InfRanges) = OneToInf() .* (first(r) .+ r) .÷ 2
cumsum(r::InfRanges) = RangeCumsum(r)
diff(r::InfRanges) = Fill(step(r),∞)
diff(r::OneToInf{T}) where T = Ones{T}(∞)
Base.@propagate_inbounds getindex(c::RangeCumsum, kr::OneToInf) = RangeCumsum(c.range[kr])
getindex(c::RangeCumsum{<:Any,<:OneToInf}, k::Integer) = k * (k+1) ÷ 2

# vcat

vcat(a::Number, r::InfRanges) = Vcat(a, r)

throw_inferror() = throw(ArgumentError("vcat is undefined with a leading infinite range"))
# vcat(r::InfRanges) = r
# vcat(r::InfRanges{T}, args::InfRanges{T}...) where {T} = throw_inferror()
# vcat(r::InfRanges, args::InfRanges...) = throw_inferror()
# vcat(r::InfRanges{T}, args::AbstractRange{T}...) where {T} = throw_inferror()
# vcat(r::InfRanges, args::AbstractRange...) = throw_inferror()
# vcat(r::InfRanges{T}, args::AbstractVector{T}...) where {T} = throw_inferror()
# vcat(r::InfRanges, args::AbstractVector...) = throw_inferror()
# vcat(r::InfRanges, ::Union{Number, AbstractVector}...) = throw_inferror()
# vcat(r::AbstractRange{T}, infr::InfRanges{T}) where {T} = Vcat(r, infr)
# vcat(r::AbstractRange, infr::InfRanges) = Vcat(r, infr)
# vcat(v::AbstractVector, infr::InfRanges) = Vcat(v, infr)

vcat(r::InfRanges{<:Number}, args::AbstractVector{<:Number}...) = throw_inferror()
vcat(r::InfRanges{<:Number}) = r
vcat(r::InfRanges{<:Number}, args::AbstractRange{<:Number}...) = throw_inferror()
vcat(r::InfRanges{<:Number}, args::InfRanges{<:Number}...) = throw_inferror()
vcat(r::InfRanges{<:Number}, ::Union{Number, AbstractVector{<:Number}}...) = throw_inferror()
vcat(r::AbstractRange{T}, infr::InfRanges{T}) where {T<:Number} = Vcat(r, infr)
vcat(r::AbstractRange{<:Number}, infr::InfRanges{<:Number}) = Vcat(r, infr)
vcat(v::AbstractVector{<:Number}, infr::InfRanges{<:Number}) = Vcat(v, infr)

###
# MemoryLayout
####

MemoryLayout(::Type{<:AbstractInfUnitRange}) = LazyLayout()
@inline MemoryLayout(::Type{<:InfStepRange}) = LazyLayout()


##
# findfirst
##

# from array.jl
function _step_findfirst(p, r::InfStepRange{T,S}) where {T,S}
    first(r) <= p.x || return nothing
    d = convert(S, p.x - first(r))
    iszero(d % step(r)) || return nothing
    return d ÷ step(r) + 1
end

for op in (:isequal, :(==))
    @eval begin
        findfirst(p::Fix2{typeof($op),T}, r::InfStepRange{T,S}) where {T,S} =
            _step_findfirst(p, r)

        findfirst(p::Fix2{typeof($op),T}, r::InfStepRange{T,<:Integer}) where {T<:Integer} =
            _step_findfirst(p, r)

        findfirst(p::Fix2{typeof($op),<:Integer}, r::InfStepRange{<:Integer,<:Integer}) =
            _step_findfirst(p, r)

        function findfirst(p::Fix2{typeof($op),<:Integer}, r::AbstractInfUnitRange{<:Integer})
            first(r) <= p.x || return nothing
            p.x - first(r) + 1
        end

        findfirst(p::Fix2{typeof($op),<:Number}, r::AbstractInfUnitRange{<:Integer}) =
            isinteger(p.x) ? findfirst($op(convert(Integer, p.x)), r) : nothing
        findfirst(p::Fix2{typeof($op),<:Number}, r::InfStepRange{<:Integer,<:Integer}) =
            isinteger(p.x) ? findfirst($op(convert(V, p.x)), r) : nothing
    end
end


FillArrays._range_convert(::Type{AbstractVector{T}}, r::InfRanges) where T = convert(AbstractVector{T}, r)

if VERSION <= v"1.9"
    permutedims(D::Diagonal{<:Any,<:InfRanges}) = D
end
copy(D::Diagonal{<:Any,<:InfRanges}) = D
broadcasted(::LazyArrayStyle{2}, ::typeof(*), a::Number, D::Diagonal{<:Any,<:InfRanges}) = a*D
broadcasted(::LazyArrayStyle{2}, ::typeof(*), D::Diagonal{<:Any,<:InfRanges}, a::Number) = D*a

function LinearAlgebra.diag(D::Diagonal{<:Any,<:InfRanges}, k::Integer = 0)
    if k == 0
        return parent(D)
    else
        return Zeros{eltype(D)}(size(D,1)) # infinite vector of zeros
    end
end
