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


"""
    imin, imax = radiuslimits(r, n::Int, m::Int)
    imin, imax = radiuslimits(r, seq1, seq2)
"""
function radiuslimits(r, n::Int, m::Int)
    d = abs(m-n)
    if m >= n
        imin, imax = max.((1:n) .- r,1), min.((1:n) .+ (r+d), m)
    else
        imin, imax = max.((1:n) .- (r+d),1), min.((1:n) .+ r, m)
    end

    imin, imax
end

radiuslimits(r,seq1, seq2) = radiuslimits(r, lastlength(seq1), lastlength(seq2))
