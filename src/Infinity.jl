## NotANumber

struct NotANumber <: Number end





"""
   Infinity()

represents infinite cardinality. Note that `Infinity <: Integer` to support
being treated as an index.
"""
struct Infinity <: Integer end

const ∞ = Infinity()

show(io::IO, y::Infinity) = print(io, "∞")

promote_rule(::Type{Infinity}, ::Type{II}) where II<:Integer = Union{Infinity,II}

convert(::Type{Float64}, ::Infinity) = Inf64
convert(::Type{Float32}, ::Infinity) = Inf32
convert(::Type{Float16}, ::Infinity) = Inf16
convert(::Type{AF}, ::Infinity) where AF<:AbstractFloat = convert(AF, Inf)


sign(y::Infinity) = 1
angle(x::Infinity) = 0

==(x::Infinity, y::Infinity) = true

isinf(::Infinity) = true
isfinite(::Infinity) = false

==(x::Infinity, y::Number) = isinf(y) && angle(y) == angle(x)
==(y::Number, x::Infinity) = x == y



isless(x::Infinity, y::Infinity) = false
isless(x::Real, y::Infinity) = isfinite(x) || sign(y) == -1
isless(x::AbstractFloat, y::Infinity) = isless(x, convert(typeof(x), y))
isless(x::Infinity, y::AbstractFloat) = false
isless(x::Infinity, y::Real) = false

+(::Infinity, ::Infinity) = ∞
+(::Number, y::Infinity) = ∞
+(::Infinity, ::Number) = ∞
-(::Infinity, ::Number) = ∞
-(x::Number, ::Infinity) = x + (-∞)

+(::Integer, y::Infinity) = ∞
+(::Infinity, ::Integer) = ∞
-(::Infinity, ::Integer) = ∞
-(x::Integer, ::Infinity) = x + (-∞)
+(::Complex, y::Infinity) = ∞
+(::Infinity, ::Complex) = ∞
-(::Infinity, ::Complex) = ∞
-(x::Complex, ::Infinity) = x + (-∞)

# ⊻ is xor
*(::Infinity) = ∞
*(::Infinity, ::Infinity) = ∞

for OP in (:fld,:cld,:div)
  @eval $OP(::Infinity, ::Number) = ∞
end


min(::Infinity, ::Infinity) = ∞
max(::Infinity, ::Infinity) = ∞
min(x::Real, ::Infinity) = x
max(::Real, ::Infinity) = ∞
min(::Infinity, x::Real) = x
max(::Infinity, ::Real) = x

≤(::Infinity, ::Infinity) = true
<(::Infinity, ::Infinity) = false
≥(::Infinity, ::Infinity) = true
>(::Infinity, ::Infinity) = false

for OP in (:<, :≤)
    @eval begin
        $OP(::Real, ::Infinity) = true
        $OP(::Infinity, ::Real) = false
    end
end


for OP in (:>, :≥)
    @eval begin
        $OP(::Real, ::Infinity) = false
        $OP(::Infinity, ::Real) = true
    end
end


# angle is π*a where a is (false==0) and (true==1)
struct OrientedInfinity{T<:Real} <: Number
    angle::T
end

OrientedInfinity{T}() where T = OrientedInfinity(zero(T))
OrientedInfinity() = OrientedInfinity{Bool}()
OrientedInfinity{T}(::Infinity) where T<:Real = OrientedInfinity{T}()
OrientedInfinity(::Infinity) = OrientedInfinity()
-(::Infinity) = OrientedInfinity(true)
+(::Infinity) = OrientedInfinity(false)


isinf(::OrientedInfinity) = true
isfinite(::OrientedInfinity) = false


promote_rule(::Type{Infinity}, ::Type{OrientedInfinity{T}}) where T = OrientedInfinity{T}
convert(::Type{OrientedInfinity{T}}, ::Infinity) where T = OrientedInfinity{T}()
convert(::Type{OrientedInfinity}, ::Infinity) = OrientedInfinity()


sign(y::OrientedInfinity{<:Integer}) = mod(y.angle,2) == 0 ? 1 : -1
angle(x::OrientedInfinity) = π*x.angle
mod(::OrientedInfinity{<:Integer}, ::Integer) = NotANumber()

function show(io::IO, y::OrientedInfinity{B}) where B<:Integer
    if sign(y) == 1
        print(io, "+∞")
    else
        print(io, "-∞")
    end
end

show(io::IO, x::OrientedInfinity) = print(io, "$(exp(im*π*x.angle))∞")

==(x::OrientedInfinity, y::Infinity) = x.angle == 0
==(y::Infinity, x::OrientedInfinity) = x.angle == 0
==(x::OrientedInfinity, y::OrientedInfinity) = x.angle == y.angle

==(x::OrientedInfinity, y::Number) = isinf(y) && angle(y) == angle(x)
==(y::Number, x::OrientedInfinity) = x == y

isless(x::OrientedInfinity{Bool}, y::OrientedInfinity{Bool}) = x.angle && !y.angle
isless(x::Number, y::OrientedInfinity{Bool}) = !y.angle && x ≠ ∞
isless(x::OrientedInfinity{Bool}, y::Number) = x.angle && y ≠ -∞

-(y::OrientedInfinity{B}) where B<:Integer = sign(y) == 1 ? OrientedInfinity(one(B)) : OrientedInfinity(zero(B))

function +(x::OrientedInfinity, y::OrientedInfinity)
    x == y || throw(ArgumentError("Angles must be the same to add ∞"))
    promote_type(typeof(x),typeof(y))(x.angle)
end

+(x::OrientedInfinity, y::Infinity) = x+OrientedInfinity(y)
+(x::Infinity, y::OrientedInfinity) = OrientedInfinity(x)+y
+(::Number, y::OrientedInfinity) = y
+(y::OrientedInfinity, ::Number) = y
-(y::OrientedInfinity, ::Number) = y
-(::Number, y::OrientedInfinity) = -y


# ⊻ is xor
*(a::OrientedInfinity{Bool}, b::OrientedInfinity{Bool}) = OrientedInfinity(a.angle ⊻ b.angle)
*(a::OrientedInfinity, b::OrientedInfinity) = OrientedInfinity(a.angle + b.angle)


*(a::Real, y::OrientedInfinity) = a > 0 ? y : (-y)
*(y::OrientedInfinity, a::Real) = a*y

*(a::Number, y::OrientedInfinity) = OrientedInfinity(y.angle+angle(a)/π)
*(y::OrientedInfinity, a::Number) = a*y

*(a::Number,y::Infinity) = a*OrientedInfinity(y)
*(y::Infinity, a::Number) = OrientedInfinity(y)*a

*(a::Integer,y::Infinity) = a*OrientedInfinity(y)
*(a::Complex,y::Infinity) = a*OrientedInfinity(y)
*(y::Infinity, a::Complex) = OrientedInfinity(y)*a
*(y::Infinity, a::Integer) = OrientedInfinity(y)*a

for OP in (:fld,:cld,:div)
  @eval $OP(y::OrientedInfinity, a::Number) = y*(1/sign(a))
end

min(x::OrientedInfinity{B}, y::OrientedInfinity{B}) where B<:Integer = sign(x) == -1 ? x : y
max(x::OrientedInfinity{B}, ::OrientedInfinity{B}) where B<:Integer = sign(x) == 1 ? x : y
min(x::Real, y::OrientedInfinity{B}) where B<:Integer = sign(y) == 1 ? x : y
min(x::OrientedInfinity{B}, y::Real) where B<:Integer = min(y,x)
max(x::Real, y::OrientedInfinity{B}) where B<:Integer = sign(y) == 1 ? y : x
max(x::OrientedInfinity{B}, y::Real) where B<:Integer = max(y,x)

for OP in (:<,:≤)
    @eval begin
        $OP(x::Real, y::OrientedInfinity{B}) where B<:Integer = sign(y) ==  1
        $OP(y::OrientedInfinity{B}, x::Real) where B<:Integer = sign(y) == -1
    end
end

for OP in (:>, :≥)
    @eval begin
        $OP(x::Real, y::OrientedInfinity{B}) where B<:Integer = sign(y) == -1
        $OP(y::OrientedInfinity{B}, x::Real) where B<:Integer = sign(y) == 1
    end
end
