module InfiniteArraysBlockArraysExt

sub_materialize(_, V, ::Tuple{BlockedOneTo{Int,<:InfRanges}}) = V
sub_materialize(::AbstractBlockLayout, V, ::Tuple{BlockedOneTo{Int,<:InfRanges}}) = V
function sub_materialize(::PaddedColumns, v::AbstractVector{T}, ax::Tuple{BlockedOneTo{Int,<:InfRanges}}) where T
    dat = paddeddata(v)
    BlockedVector(Vcat(sub_materialize(dat), Zeros{T}(∞)), ax)
end

end # module