for typ in (:Fill, :Zeros, :Ones)
    @eval begin
        Base.IndexStyle(::Type{<:$typ{<:Any,2,Tuple{Infinity,Infinity}}}) = Base.IndexCartesian()
        Base.IndexStyle(::Type{<:$typ{<:Any,2,Tuple{Infinity,Int}}}) = Base.IndexCartesian()
    end
end
# axes(::Fill{<:Any,1,Tuple{Infinity}}) = tuple(OneToInf())
