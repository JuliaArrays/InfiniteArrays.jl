##
# hash
#
# `Base.hash(::AbstractArray, ::UInt)` hashes the endpoints of the axes and then walks the
# entries backwards from the last index, which neither terminates nor is well-defined for an
# infinite range. We follow the same recipe, but take the entries from a finite window at the
# start instead. Equal ranges have equal axes and equal entries, so this preserves the contract
# `isequal(a,b) ⟹ hash(a) == hash(b)`; the price is that ranges which first differ outside the
# window collide.
#
# Only the ranges defined here are covered. The other infinite arrays of the ecosystem are hashed
# by the packages that own them, using this same scheme — keep the seeds and the window size in
# sync with `FillArrays`, `ArrayLayouts` and `BlockArrays` so that arrays which compare equal
# across packages keep hashing alike.
##

const hash_infinity_seed = UInt === UInt64 ? 0x3a4c1b6d9e2f5087 : 0x9e2f5087
const hash_infarray_seed = UInt === UInt64 ? 0x5d7e3f11a8c6b024 : 0xa8c6b024

# number of entries hashed along each dimension of an infinite array
const hash_nentries = 8

# The endpoint of an infinite axis is an infinity, which `Base` hashes only from Infinities
# v0.1.13 on, so hash the direction it points in instead: infinities pointing the same way compare
# equal and hence have to hash alike.
_hash_endpoint(x, h::UInt) = isinf(x) ? hash(signbit(x), h + hash_infinity_seed) : hash(x, h)

"""
    _hash_indices(ax)

The indices along the axis `ax` at which entries get hashed: the first `hash_nentries` of them, or
all of `ax` when it is shorter than that. A bi-infinite axis has no first index, so there the
window starts at zero instead.
"""
function _hash_indices(ax::AbstractUnitRange)
    a = first(ax)
    isinf(a) && return 0:hash_nentries-1
    b = a + hash_nentries - 1
    a:(isinf(last(ax)) ? b : min(last(ax), b))
end
_hash_indices(ax) = Iterators.take(ax, hash_nentries)

function _hash_infarray(A::AbstractArray, h::UInt)
    h += hash_infarray_seed
    # the axes are arrays in their own right, so hash their endpoints instead of the axes
    for ax in axes(A)
        h = _hash_endpoint(first(ax), h)
    end
    for ax in axes(A)
        h = _hash_endpoint(last(ax), h)
    end
    for I in Iterators.product(map(_hash_indices, axes(A))...)
        h = hash(A[I...], h)
    end
    h
end

# `length` is the product of the axis lengths, and `0 * ℵ₀` throws, so ask the axes one at a time
_hasinfaxes(A::AbstractArray) = any(ax -> isinf(length(ax)), axes(A))

"""
    InfRangeLike

The infinite ranges defined here together with the `Base` wrappers used to turn them into axes:
`Slice(OneToInf())` and `IdentityUnitRange(1:∞)` compare equal to the range they wrap and so have
to hash alike, and no other package is in a position to define that.
"""
const InfRangeLike = Union{InfRanges,
                           Slice{<:AbstractInfUnitRange},
                           IdentityUnitRange{<:AbstractInfUnitRange}}

# the `LinearAlgebra` matrices that are infinite when built out of an infinite range
const InfRangeWrapper = Union{Diagonal{<:Any,<:InfRangeLike},
                              Bidiagonal{<:Any,<:InfRangeLike},
                              Tridiagonal{<:Any,<:InfRangeLike},
                              SymTridiagonal{<:Any,<:InfRangeLike}}

"""
    InfRangeWrappers

The infinite ranges together with the `LinearAlgebra` and `Base` wrappers of them. Such a wrapper
can only be covered here: `FillArrays` and `ArrayLayouts` cover the wrappers of their own arrays,
but neither is in a position to recognise an infinite range. Not every instance is infinite —
`view(1:∞, 1:5)` is not — so `hash` checks at run time and leaves the finite ones to `Base`.
"""
const InfRangeWrappers = Union{InfRangeLike,
                               InfRangeWrapper,
                               HermOrSym{<:Any,<:InfRangeWrapper},
                               UpperOrLowerTriangular{<:Any,<:InfRangeWrapper},
                               AdjOrTrans{<:Any,<:Union{InfRangeLike,InfRangeWrapper}},
                               SubArray{<:Any,<:Any,<:Union{InfRangeLike,InfRangeWrapper}}}

hash(A::InfRangeWrappers, h::UInt) =
    _hasinfaxes(A) ? _hash_infarray(A, h) : invoke(hash, Tuple{AbstractArray,UInt}, A, h)

# `Base` compares ranges structurally in `==` but defines no `isequal` for them, so `isequal` falls
# back to walking the entries and never terminates. Without this, infinite ranges cannot be used as
# `Dict` or `Set` keys even once they hash. Like `Base.isequal(::AbstractArray, ::AbstractArray)` we
# compare the axes first: an offset range such as `IdentityUnitRange(3:∞)` holds the same entries as
# `3:∞` but is not `isequal` to it, and hashes differently.
isequal(r::InfRangeLike, s::InfRangeLike) =
    axes(r) == axes(s) && isequal(first(r), first(s)) && isequal(step(r), step(s)) &&
        isequal(last(r), last(s))
