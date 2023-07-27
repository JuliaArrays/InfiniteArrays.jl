module InfiniteArraysStatisticsExt

using InfiniteArrays
using InfiniteArrays: InfRanges
using Statistics

Statistics.mean(r::InfRanges{<:Real}) = last(r)
Statistics.median(r::InfRanges{<:Real}) = last(r)

end
