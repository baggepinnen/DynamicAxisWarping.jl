Base.@propagate_inbounds Base.getindex(v::AbstractVector, ::typeof(!), i::Int) = v[i]
Base.@propagate_inbounds Base.getindex(v::AbstractVector, ::typeof(!), i::Int) = @view v[:,i]

@inline function indmin3(a,b,c,i,j)
    if a <= b
        if a <= c
            return 1,i-1,j-1
        else
            return 3,i,j-1
        end
    else
        if b <= c
            return 2,i-1,j
        else
            return 3,i,j-1
        end
    end
end
