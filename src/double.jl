
# generic interface
export shift

# interface for compact sequences
export zero_padding, CompactSequence
export support

# interface for periodic extensions
export PeriodicExtension
export period


"""
A `DoublyInfiniteVector` is a sequence of infinite length that can be indexed
with any integer, including both positive and negative ones.

It is also referred to as a bi-infinite sequence or two-way infinite sequence.
"""
abstract type DoublyInfiniteVector{T} <: InfVector{T}
end


getindex(s::DoublyInfiniteVector, r::Range) = [s[i] for i in r]

# The parahermitian conjugate of a sequence is the conjugate of its reverse
parahermitian_conjugate(s::DoublyInfiniteVector) = conj(reverse(s))

adjoint(s::DoublyInfiniteVector) = conj(s)


"The supertype of doubly infinite sequences obtained by extending finite data."
abstract type ExtensionSequence{T} <: DoublyInfiniteVector{T}
end

"The length of the data vector that is stored in the sequence."
data_length(s::ExtensionSequence) = length(data(s))

"The data_left index of the embedded data in the sequence."
data_leftindex(s::ExtensionSequence) = offset(s)
"The data_right index of the embedded data in the sequence."
data_rightindex(s::ExtensionSequence) = offset(s) + data_length(s) - 1
"The support of the embedded data in the sequence."
data_support(s::ExtensionSequence) = data_leftindex(s):data_rightindex(s)

"The convex hull of the supports of the embedded data in two sequences."
function data_support(s1::ExtensionSequence, s2::ExtensionSequence)
    l = min(data_leftindex(s1), data_leftindex(s2))
    r = max(data_rightindex(s1), data_rightindex(s2))
    l:r
end

"Shift the sequence `k` positions to the right."
shift(s::ExtensionSequence, k::Int) = similar(s, data(s), offset(s)+k)


# We can do arithmetics generically by performing operations on the underlying
# data. Concrete sequences should implement `align` which aligns the data so
# that this is possible.
for op in (:+, :-)
    @eval function $op(s1::S, s2::S) where S <: ExtensionSequence
        a, b = align(s1, s2)
        similar(s1, $op(data(a), data(b)), offset(a))
    end
end

for op in (:+, :-, :*, :/)
    @eval ($op)(s::ExtensionSequence, x::Number) = similar(s, ($op)(data(s),x), offset(s))
end

for op in (:+, :-, :*)
    @eval ($op)(x::Number, s::ExtensionSequence) = ($op)(s, x)
end

-(s::ExtensionSequence) = similar(s, -data(s), offset(s))



################
## Sequences with compact support, or extension by zero padding
################

# For documentation see constructor below
struct CompactSequence{T} <: ExtensionSequence{T}
    # We store the non-zero coefficients as a regular vector
    v       ::  Vector{T}
    # offset stores the index at which the data starts, default is 1
    offset  ::  Int
end


"""
    CompactSequence(v[, offset])

A `CompactSequence` is a compactly supported sequence, i.e., a sequence with
a finite number of nonzero elements. Its support is an interval ``[i,j]`` where
both `i` and `j` are integers and `j >= i`. The sequence elements are zero
outside that interval.

The `CompactSequence` stores a vector of length `j-i+1`. It can be thought of
as an extension of this vector to doubly infinite sequences by zero padding.
The optional `offset` in the constructor is the left endpoint `i` of the support
of the sequence. By default, it equals the first valid index of `v`.
"""
CompactSequence(v::AbstractVector) = CompactSequence(v, first(eachindex(v)))

"Extend the given vector to a doubly infinite sequency by zero padding."
zero_padding(v::AbstractVector, optional...) = CompactSequence(v, optional...)

data(s::CompactSequence) = s.v
offset(s::CompactSequence) = s.offset

similar(s::CompactSequence, data::AbstractVector, optional...) = CompactSequence(data, optional...)

"A range of indices that includes all non-zero elements of the sequence."
support(s::CompactSequence) = data_support(s)

# For internal use:
# - map the index k of a sequence into an index l of the data
mapindex(s::CompactSequence, k) = k - offset(s) + 1
# - and vice-versa
imapindex(s::CompactSequence, l) = l + offset(s) - 1


# We override getindex to return zero outside our embedded vector.
getindex(s::CompactSequence, k::Int) =
    k < data_leftindex(s) || k > data_rightindex(s) ? zero(eltype(s)) : getindex(data(s), mapindex(s, k))

# Reverse the sequence in time
reverse(s::CompactSequence) = CompactSequence(flipdim(data(s), 1), -data_rightindex(s))

# Take element-wise conjugates
conj(s::CompactSequence) = CompactSequence(conj(data(s)), data_leftindex(s))

Base.widen(s::CompactSequence) = CompactSequence(widen(data(s)), offset(s))

# Return two sequences that have the same support, by zero padding
function align(s1::CompactSequence, s2::CompactSequence)
    if support(s1) == support(s2)
        # they already align, do nothing
        s1, s2
    else
        #they do not align, make new sequences
        supp = data_support(s1, s2)
        CompactSequence(s1[supp], first(supp)), CompactSequence(s2[supp], first(supp))
    end
end


################
## Periodic extension
################

# For documentation see constructor below
struct PeriodicExtension{T} <: ExtensionSequence{T}
    v       ::  Vector{T}
    offset  ::  Int
end

"""
    PeriodicExtension(v[, offset])

Extend a given vector `v` using periodicity to a doubly infinite sequence.

The `offset` is optional and by default equals the first index of `v`. The sequence
is defined by
``
h[i] =
``
"""
PeriodicExtension(v::AbstractVector) = PeriodicExtension(v, first(eachindex(v)))

data(s::PeriodicExtension) = s.v
offset(s::PeriodicExtension) = s.offset

period(s::PeriodicExtension) = data_length(s)

similar(s::PeriodicExtension, data::AbstractVector, optional...) =
    PeriodicExtension(data, optional...)

# For internal use:
# - map the index k of a sequence into an index l of the data
mapindex(s::PeriodicExtension, k) = mod(k - offset(s), period(s)) + 1
# - and vice-versa
imapindex(s::PeriodicExtension, l) = l + offset(s) - 1

getindex(s::PeriodicExtension, k::Int) = getindex(data(s), mapindex(s, k))

# Return two periodic sequences that are aligned (have the same length and offset)
function align(s1::PeriodicExtension, s2::PeriodicExtension)
    @assert period(s1) == period(s2)
    if offset(s1) == offset(s2)
        s1, s2
    else
        s1, PeriodicExtension(circshift(data(s2), offset(s2)-offset(s1)), offset(s1))
    end
end

################
## Symmetric extensions
################
