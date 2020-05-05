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


function radiuslimits(r,n, m)
    max.((1:n) .- r,1), min.((1:n) .+ r, m)
end
