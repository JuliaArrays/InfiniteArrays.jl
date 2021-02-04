## NotANumber

struct NotANumber <: Number end





"""
   Infinity()

represents infinite cardinality. Note that `Infinity <: Integer` to support
being treated as an index.
"""
struct Infinity <: Integer end

const ∞ = Infinity()

show(io::IO, ::Infinity) = print(io, "∞")
string(::Infinity) = "∞"

promote_rule(::Type{Infinity}, ::Type{II}) where II<:Integer = Union{Infinity,II}

convert(::Type{Float64}, ::Infinity) = Inf64
convert(::Type{Float32}, ::Infinity) = Inf32
convert(::Type{Float16}, ::Infinity) = Inf16
convert(::Type{AF}, ::Infinity) where AF<:AbstractFloat = convert(AF, Inf)


sign(y::Infinity) = 1
angle(x::Infinity) = 0

==(x::Infinity, y::Infinity) = true

one(::Type{Infinity}) = 1
zero(::Infinity) = 0

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

-(::Infinity, ::Infinity) = NotANumber()

# ⊻ is xor
*(::Infinity) = ∞
*(::Infinity, ::Infinity) = ∞



for OP in (:fld,:cld,:div)
  @eval begin
    $OP(::Infinity, ::Real) = ∞
    $OP(::Infinity, ::Infinity) = NotANumber()
  end
end

div(::T, ::Infinity) where T<:Real = zero(T)
fld(x::T, ::Infinity) where T<:Real = signbit(x) ? -one(T) : zero(T)
cld(x::T, ::Infinity) where T<:Real = signbit(x) ? zero(T) : one(T)

mod(::Infinity, ::Infinity) = NotANumber()
mod(::Infinity, ::Real) = NotANumber()
function mod(x::Real, ::Infinity) 
    x ≥ 0 || throw(ArgumentError("mod(x,∞) is unbounded for x < 0"))
    x
end

min(::Infinity, ::Infinity) = ∞
max(::Infinity, ::Infinity) = ∞
min(x::Real, ::Infinity) = x
max(::Real, ::Infinity) = ∞
min(::Infinity, x::Real) = x
max(::Infinity, ::Real) = ∞

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


struct SignedInfinity <: Integer
    signbit::Bool
end

SignedInfinity() = SignedInfinity(false)
SignedInfinity(::Infinity) = SignedInfinity()
SignedInfinity(x::SignedInfinity) = x

-(::Infinity) = SignedInfinity(true)
+(::Infinity) = ∞

isinf(::SignedInfinity) = true
isfinite(::SignedInfinity) = false

promote_rule(::Type{Infinity}, ::Type{SignedInfinity}) = SignedInfinity
convert(::Type{SignedInfinity}, ::Infinity) = SignedInfinity(false)

signbit(y::SignedInfinity) = y.signbit
sign(y::SignedInfinity) = 1-2signbit(y)
angle(x::SignedInfinity) = π*signbit(x)
mod(::SignedInfinity, ::SignedInfinity) = NotANumber()
mod(::SignedInfinity, ::Real) = NotANumber()
function mod(x::Real, y::SignedInfinity) 
    signbit(x) == signbit(y) || throw(ArgumentError("mod($x,$y) is unbounded"))
    x
end

string(y::SignedInfinity) = signbit(y) ? "-∞" : "+∞"
show(io::IO, y::SignedInfinity) = print(io, string(y))

==(x::SignedInfinity, y::Infinity) = !x.signbit
==(y::Infinity, x::SignedInfinity) = !x.signbit
==(x::SignedInfinity, y::SignedInfinity) = x.signbit == y.signbit

==(x::SignedInfinity, y::Number) = isinf(y) && signbit(y) == signbit(x)
==(y::Number, x::SignedInfinity) = x == y

isless(x::SignedInfinity, y::SignedInfinity) = signbit(x) && !signbit(y)
for Typ in (:Number, :Real, :Integer, :AbstractFloat)
    @eval begin
        isless(x::SignedInfinity, y::$Typ) = signbit(x) && y ≠ -∞
        isless(x::$Typ, y::SignedInfinity) = !signbit(y) && x ≠ ∞
        +(::$Typ, y::SignedInfinity) = y
        +(y::SignedInfinity, ::$Typ) = y
        -(y::SignedInfinity, ::$Typ) = y
        -(::$Typ, y::SignedInfinity) = -y
        function *(a::$Typ, y::SignedInfinity) 
            iszero(a) && throw(ArgumentError("Cannot multiply $a * $y"))
            a > 0 ? y : (-y)
        end
    end
end

≤(::SignedInfinity, ::Infinity) = true
≤(::Infinity, s::SignedInfinity) = !signbit(s)
<(s::SignedInfinity, ::Infinity) = signbit(s)
<(::Infinity, ::SignedInfinity) = false
≥(s::SignedInfinity, ::Infinity) = !signbit(s)
≥(::Infinity, ::SignedInfinity) = true
>(::SignedInfinity, ::Infinity) = false
>(::Infinity, s::SignedInfinity) = signbit(s)




function -(::Infinity, y::SignedInfinity) 
    signbit(y) || throw(ArgumentError("Cannot subtract ∞ from ∞"))
    ∞
end

function -(x::SignedInfinity, ::Infinity) 
    signbit(x) || throw(ArgumentError("Cannot subtract ∞ from ∞"))
    x
end

function -(x::SignedInfinity, y::SignedInfinity) 
    signbit(x) == !signbit(y) || throw(ArgumentError("Cannot subtract ∞ from ∞"))
    x
end

-(y::SignedInfinity) = SignedInfinity(!y.signbit)

function +(x::SignedInfinity, y::SignedInfinity)
    x == y || throw(ArgumentError("Angles must be the same to add ∞"))
    x
end

+(x::SignedInfinity, y::Infinity) = x+SignedInfinity(y)
+(x::Infinity, y::SignedInfinity) = SignedInfinity(x)+y

# ⊻ is xor
*(a::SignedInfinity, b::SignedInfinity) = SignedInfinity(signbit(a) ⊻ signbit(b))
*(a::Infinity, b::SignedInfinity) = SignedInfinity(a)*b
*(a::SignedInfinity, b::Infinity) = a*SignedInfinity(b)

*(a::Integer, y::Infinity) = a*SignedInfinity(y)
*(y::Infinity, a::Integer) = SignedInfinity(y)*a

*(a::Real, y::Infinity) = a*SignedInfinity(y)
*(y::Infinity, a::Real) = SignedInfinity(y)*a

*(y::SignedInfinity, a::Real) = a*y
*(y::SignedInfinity, a::Integer) = a*y

<(x::SignedInfinity, y::SignedInfinity) = signbit(x) & !signbit(y)
≤(x::SignedInfinity, y::SignedInfinity) = signbit(x) | !signbit(y)

for OP in (:<,:≤)
    @eval begin
        $OP(x::Real, y::SignedInfinity) = !signbit(y)
        $OP(y::SignedInfinity, x::Real) = signbit(y)
    end
end


min(x::SignedInfinity, y::SignedInfinity) = SignedInfinity(x.signbit | y.signbit)
max(x::SignedInfinity, y::SignedInfinity) = SignedInfinity(x.signbit & y.signbit)
min(x::Real, y::SignedInfinity) = y.signbit ? x : y
max(x::Real, y::SignedInfinity) = y.signbit ? y : x
min(x::SignedInfinity, y::Real) = x.signbit ? x : y
max(x::SignedInfinity, y::Real) = x.signbit ? y : x
min(x::SignedInfinity, ::Infinity) = x
max(::SignedInfinity, ::Infinity) = ∞
min(::Infinity, x::SignedInfinity) = x
max(::Infinity, x::SignedInfinity) = ∞



######
# OrientedInfinity
#######

# angle is π*a where a is (false==0) and (true==1)
struct OrientedInfinity{T<:Real} <: Number
    angle::T
end

OrientedInfinity{T}() where T = OrientedInfinity(zero(T))
OrientedInfinity() = OrientedInfinity{Bool}()
OrientedInfinity{T}(::Infinity) where T<:Real = OrientedInfinity{T}()
OrientedInfinity(::Infinity) = OrientedInfinity()
OrientedInfinity{T}(x::SignedInfinity) where T<:Real = OrientedInfinity{T}(signbit(x))
OrientedInfinity(x::SignedInfinity) = OrientedInfinity(signbit(x))



isinf(::OrientedInfinity) = true
isfinite(::OrientedInfinity) = false


promote_rule(::Type{Infinity}, ::Type{OrientedInfinity{T}}) where T = OrientedInfinity{T}
convert(::Type{OrientedInfinity{T}}, ::Infinity) where T = OrientedInfinity{T}()
convert(::Type{OrientedInfinity}, ::Infinity) = OrientedInfinity()
convert(::Type{OrientedInfinity{T}}, x::SignedInfinity) where T = OrientedInfinity{T}(x)
convert(::Type{OrientedInfinity}, x::SignedInfinity) = OrientedInfinity(x)


sign(y::OrientedInfinity{<:Integer}) = mod(y.angle,2) == 0 ? 1 : -1
angle(x::OrientedInfinity) = π*x.angle
mod(::OrientedInfinity{<:Integer}, ::Integer) = NotANumber()


show(io::IO, x::OrientedInfinity) = print(io, "$(exp(im*π*x.angle))∞")

==(x::OrientedInfinity, y::Infinity) = x.angle == 0
==(y::Infinity, x::OrientedInfinity) = x.angle == 0
==(x::OrientedInfinity, y::SignedInfinity) = x.angle == signbit(y)
==(y::SignedInfinity, x::OrientedInfinity) = x.angle == signbit(y)
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
+(x::OrientedInfinity, y::SignedInfinity) = x+OrientedInfinity(y)
+(x::SignedInfinity, y::OrientedInfinity) = OrientedInfinity(x)+y
+(::Number, y::OrientedInfinity) = y
+(y::OrientedInfinity, ::Number) = y
-(y::OrientedInfinity, ::Number) = y
-(::Number, y::OrientedInfinity) = -y


# ⊻ is xor
*(a::OrientedInfinity{Bool}, b::OrientedInfinity{Bool}) = OrientedInfinity(a.angle ⊻ b.angle)
*(a::OrientedInfinity, b::OrientedInfinity) = OrientedInfinity(a.angle + b.angle)
*(a::Infinity, b::OrientedInfinity) = OrientedInfinity(a)*b
*(a::OrientedInfinity, b::Infinity) = a*OrientedInfinity(b)
*(a::SignedInfinity, b::OrientedInfinity) = OrientedInfinity(a)*b
*(a::OrientedInfinity, b::SignedInfinity) = a*OrientedInfinity(b)

*(a::Real, y::OrientedInfinity) = a > 0 ? y : (-y)
*(y::OrientedInfinity, a::Real) = a*y

*(a::Number, y::OrientedInfinity) = OrientedInfinity(y.angle+angle(a)/π)
*(y::OrientedInfinity, a::Number) = a*y

*(a::Number, y::Infinity) = a*OrientedInfinity(y)
*(y::Infinity, a::Number) = OrientedInfinity(y)*a
*(y::SignedInfinity, a::Number) = OrientedInfinity(y)*a

*(a::Complex, y::Infinity) = a*OrientedInfinity(y)
*(y::Infinity, a::Complex) = OrientedInfinity(y)*a

*(a::Complex,y::SignedInfinity) = a*OrientedInfinity(y)
*(y::SignedInfinity, a::Complex) = OrientedInfinity(y)*a

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

##
# Checked
##

Base.Checked.checked_sub(::Integer, ::Infinity) = -∞
Base.Checked.checked_sub(::Infinity, ::Integer) = ∞
Base.Checked.checked_add(::Integer, ::Infinity) = ∞
Base.Checked.checked_add(::Infinity, ::Integer) = ∞

Base.Checked.checked_sub(::Integer, x::SignedInfinity) = -x
Base.Checked.checked_sub(x::SignedInfinity, ::Integer) = x
Base.Checked.checked_add(::Integer, x::SignedInfinity) = x
Base.Checked.checked_add(x::SignedInfinity, ::Integer) = x


Base.to_index(::Infinity) = ∞


Base.hash(::Infinity) = 0x020113134b21797f # made up