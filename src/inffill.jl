for typ in (:Fill, :Zeros, :Ones)
    @eval begin
        Base.IndexStyle(::Type{<:$typ{<:Any,2,Tuple{Infinity,Infinity}}}) = Base.IndexCartesian()
        Base.IndexStyle(::Type{<:$typ{<:Any,2,Tuple{Infinity,Int}}}) = Base.IndexCartesian()
    end
end
# axes(::Fill{<:Any,1,Tuple{Infinity}}) = tuple(OneToInf())

# Lazy Broadacasting
for typ in (:Ones, :Zeros, :Fill)
    @eval begin
        BroadcastStyle(::Type{$typ{T,N,NTuple{N,Infinity}}}) where {T,N} = LazyArrayStyle{N}()
        BroadcastStyle(::Type{$typ{T,2,Tuple{Int,Infinity}}}) where {T} = LazyArrayStyle{2}()
        BroadcastStyle(::Type{$typ{T,2,Tuple{Infinity,Int}}}) where {T} = LazyArrayStyle{2}()
    end
end

BroadcastStyle(::Type{Eye{T,NTuple{2,Infinity}}}) where {T} = LazyArrayStyle{2}()
BroadcastStyle(::Type{Eye{T,Tuple{Int,Infinity}}}) where {T} = LazyArrayStyle{2}()
BroadcastStyle(::Type{Eye{T,Tuple{Infinity,Int}}}) where {T} = LazyArrayStyle{2}()
