
export zero_padding, CompactSequence
export leftmost, rightmost

"""
A `DoublyInfiniteVector` is a sequence of infinite length that can be indexed
with any integer, including both positive and negative ones.

It is also referred to as a bi-infinite sequence or two-way infinite sequence.
"""
abstract type DoublyInfiniteVector{T} <: InfVector{T}
end


################
## Sequences with compact support, or extension by zero padding
################

# For documentation see constructor below
struct CompactSequence{T} <: DoublyInfiniteVector{T}
    # We store the non-zero coefficients as a regular vector
    data    ::  Vector{T}
    # offset stores the index at which the data starts, default is 1
    offset  ::  Int
end


"""
    CompactSequence(data[, offset])

A `CompactSequence` is a compactly supported sequence, i.e., a sequence with
a finite number of nonzero elements. Its support is an interval ``[i,j]`` where
both `i` and `j` are integers and `j >= i`. The sequence elements are zero
outside that interval.

The `CompactSequence` stores a vector of length `j-i+1`. It can be thought of
as an extension of this vector to doubly infinite sequences by zero padding.
The optional `offset` in the constructor is the left endpoint `i` of the support
of the sequence. By default, it equals the first valid index of the given vector.
"""
CompactSequence(data::AbstractVector) = CompactSequence(data, first(eachindex(data)))

"Extend the given vector to a doubly infinite sequency by zero padding."
zero_padding(data::AbstractVector, optional...) = CompactSequence(data, optional...)

data(s::CompactSequence) = s.data
offset(s::CompactSequence) = s.offset

"Return the leftmost index of the support of the sequence."
leftmost(s::CompactSequence) = s.offset
"Return the rightmost index of the support of the sequence."
rightmost(s::CompactSequence) = s.offset + length(s.data) - 1

"The length of the vector that is stored in the compact sequence."
datalength(s::CompactSequence) = length(data(s))

"A range of indices that includes all non-zero elements of the sequence."
support(s::CompactSequence) = leftmost(s):rightmost(s)

"Shift a compact sequence `k` positions to the right."
shift(s::CompactSequence, k::Int) = CompactSequence(data(s), offset(s)+k)

# For internal use:
# - map the index k of a sequence into an index l of the data
mapindex(s::CompactSequence, k) = k - offset(s) + 1
# - and vice-versa
imapindex(s::CompactSequence, l) = l + offset(s) - 1


# We override getindex to return zero outside our embedded vector.
getindex(s::CompactSequence, k::Int) =
    k < leftmost(s) || k > rightmost(s) ? zero(eltype(s)) : getindex(data(s), mapindex(s, k))

getindex(s::CompactSequence, r::Range) = [s[i] for i in r]

# Define in the future: infvector from indexing a doubly infinite sequence with an infinite range
#getindex(s::CompactSequence, r::InfUnitRange)

# Reverse the sequence in time
reverse(s::CompactSequence) = CompactSequence(flipdim(data(s), 1), -rightmost(s))

# Take element-wise conjugates
conj(s::CompactSequence) = CompactSequence(conj(data(s)), leftmost(s))

# We define the parahermitian conjugate of a sequence as the conjugate of its reverse
parahermitian_conjugate(s::CompactSequence) = conj(reverse(s))

Base.widen(s::CompactSequence) = CompactSequence(widen(data(s)), offset(s))

for op in (:+, :-, :*)
    @eval ($op)(x::Number, s::CompactSequence) = CompactSequence(($op)(s.a,x),s.offset)
    @eval ($op)(s::CompactSequence, x::Number) = ($op)(x, s)
end

for op in (:+, :-)
    @eval function $op(s1::CompactSequence, s2::CompactSequence)
        l = min(leftmost(s1), leftmost(s2))
        r = max(rightmost(s1), rightmost(s2))
        # This is not the most efficient method in all cases, but it is generic
        CompactSequence([s1[i] + s2[i] for i in l:r], l)
    end
end

-(s::CompactSequence) = CompactSequence(-data(s), offset(s))

/(s::CompactSequence, x::Number) = CompactSequence(data(s)/x, offset(s))


################
## Periodic extension
################
