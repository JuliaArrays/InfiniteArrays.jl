# This gets called when infinit number of columns
print_matrix_row(io::IO,
        X::AbstractVecOrMat, A::Vector,
        i::Integer, cols::AbstractVector{<:Infinity}, sep::AbstractString) = nothing


print_matrix_vdots(io::IO, vdots::AbstractString,
        A::Vector, sep::AbstractString, M::Integer, ::NotANumber) = nothing
