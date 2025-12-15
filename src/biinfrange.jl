"""
BiInfUnitRange()

Represent -∞:∞ with offset indexing
"""
struct BiInfUnitRange{T<:Real} <: AbstractInfUnitRange{T} end

BiInfUnitRange() = BiInfUnitRange{Int}()

AbstractArray{T}(a::BiInfUnitRange) where T<:Real = BiInfUnitRange{T}(a)
AbstractVector{T}(a::BiInfUnitRange) where T<:Real = BiInfUnitRange{T}(a)

unitrange(a::BiInfUnitRange) = a
Base.has_offset_axes(::BiInfUnitRange) = true

getindex(v::BiInfUnitRange{T}, i::Integer) where T = convert(T, i)
getindex(v::BiInfUnitRange{T}, i::RealInfinity) where T = i
axes(::BiInfUnitRange) = (BiInfUnitRange(),)
first(::BiInfUnitRange) = -∞
show(io::IO, ::BiInfUnitRange{Int}) = print(io, "BiInfUnitRange()")

getindex(r::BiInfUnitRange{T}, s::AbstractUnitRange{<:Integer}) where T = convert(AbstractVector{T}, s)



function Base._print_matrix(io, @nospecialize(X::AbstractVecOrMat), pre, sep, post, hdots, vdots, ddots, hmod, vmod, rowsA::BiInfUnitRange, colsA)
    hmod, vmod = Int(hmod)::Int, Int(vmod)::Int
    ncols, idxlast = length(colsA), last(colsA)
    if !(get(io, :limit, false)::Bool)
        screenheight = screenwidth = typemax(Int)
    else
        sz = displaysize(io)::Tuple{Int,Int}
        screenheight, screenwidth = sz[1] - 4, sz[2]
    end
    screenwidth -= length(pre)::Int + length(post)::Int
    presp = repeat(" ", length(pre)::Int)  # indent each row to match pre string
    postsp = ""
    @assert textwidth(hdots) == textwidth(ddots)
    sepsize = length(sep)::Int
    m, n = length(rowsA), length(colsA)
    # To figure out alignments, only need to look at as many rows as could
    # fit down screen. If screen has at least as many rows as A, look at A.
    # If not, then we only need to look at the first and last chunks of A,
    # each half a screen height in size.
    quarterheight = div(screenheight,4)
    halfheight = 2quarterheight+1
    rowsA = -quarterheight:quarterheight
    # Similarly for columns, only necessary to get alignments for as many
    # columns as could conceivably fit across the screen
    maxpossiblecols = div(screenwidth, 1+sepsize)
    if n > maxpossiblecols
        colsA = [colsA[(0:maxpossiblecols-1) .+ firstindex(colsA)]; colsA[(end-maxpossiblecols+1):end]]
    else
        colsA = [colsA;]
    end
    A = alignment(io, X, rowsA, colsA, screenwidth, screenwidth, sepsize, ncols)
    # Nine-slicing is accomplished using print_matrix_row repeatedly
    if n <= length(A) # rows don't fit, cols do, so only vertical ellipsis
        print_matrix_vdots(io, vdots, A, sep, vmod, 1, false)
        print(io, postsp * '\n')
        for i in rowsA
            print(io, i == first(rowsA) ? pre : presp)
            i == 0 && print(io, "\e[1m")
            print_matrix_row(io, X,A,i,colsA,sep,idxlast)
            i == 0 && print(io, "\e[0m")
            print(io, i == last(rowsA) ? post : postsp)
            if i != rowsA[end] || i == rowsA[halfheight]; println(io); end
            if i == rowsA[halfheight]
                print(io, i == first(rowsA) ? pre : presp)
                print_matrix_vdots(io, vdots, A, sep, vmod, 1, false)
                print(io, i == last(rowsA) ? post : postsp * '\n')
            end
        end
    else # neither rows nor cols fit, so use all 3 kinds of dots
        c = div(screenwidth-length(hdots)::Int+1,2)+1
        Ralign = reverse(alignment(io, X, rowsA, reverse(colsA), c, c, sepsize, ncols))
        c = screenwidth - sum(map(sum,Ralign)) - (length(Ralign)-1)*sepsize - length(hdots)::Int
        Lalign = alignment(io, X, rowsA, colsA, c, c, sepsize, ncols)
        r = mod((length(Ralign)-n+1),vmod) # where to put dots on right half
        for i in rowsA
            print(io, i == first(rowsA) ? pre : presp)
            print_matrix_row(io, X,Lalign,i,colsA[1:length(Lalign)],sep,idxlast)
            print(io, (i - first(rowsA)) % hmod == 0 ? hdots : repeat(" ", length(hdots)::Int))
            print_matrix_row(io, X,Ralign,i,(n-length(Ralign)).+colsA,sep,idxlast)
            print(io, i == last(rowsA) ? post : postsp)
            if i != rowsA[end] || i == rowsA[halfheight]; println(io); end
            if i == rowsA[halfheight]
                print(io, i == first(rowsA) ? pre : presp)
                print_matrix_vdots(io, vdots, Lalign, sep, vmod, 1, true)
                print(io, ddots)
                print_matrix_vdots(io, vdots, Ralign, sep, vmod, r, false)
                print(io, i == last(rowsA) ? post : postsp * '\n')
            end
        end
    end
    if isempty(rowsA)
        print(io, pre)
        print(io, vdots)
        length(colsA) > 1 && print(io, "    ", ddots)
        print(io, post)
    end
end
