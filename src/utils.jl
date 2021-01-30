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

"""
    inds = align_signals(s::AbstractVector{<:AbstractVector}, master = argmax(length.(s)); method = :dtw)

Compute a set of indices such that `s[i][inds[i]]` is optimally aligned to `s[master]`. 

# Arguments:
- `s`: A vector of signals to align
- `master`: Index of the signal used as reference.
- `method`: `:dtw` uses the warping paths from dtw between `s[master]` and `s[i]`. `:xcorr` uses `DSP.finddelay` which internally computes the cross correlation between signals, which often results in a slight misalignment.  
"""
function align_signals(s::AbstractVector{<:AbstractVector}, master::Integer=argmax(length.(s)); method=:dtw)
    inds = UnitRange.(eachindex.(s))
    # find delays to align with master
    d = map(s) do si
        si === s[master] && (return 0)
        if method ∈ (:xcorr, :crosscorr, :dsp)
            SlidingDistancesBase.DSP.finddelay(s[master], si) # suboptimal because xcorr does not do exactly what we want
        elseif method ∈ (:dtw, :DTW)
            d,i1,i2 = dtw(si, s[master])
            round(Int, median(i2-i1))
        else
            throw(ArgumentError("Unknown method"))
        end
    end

    # find left and right (virtual) zero padding
    lp = maximum(d)
    rp = maximum(length(s[master]) .- (length.(s) .+ d))
    
    # New window length
    wl = length(inds[master]) - lp - rp
    
    # trim individual index sets to fit into new master window
    for i in eachindex(inds)
        start = max(1, 1+lp-d[i])
        stop = min(length(s[i]),start+wl-1)
        inds[i] = start : stop
    end
    @assert all(length.(inds) .== length(inds[1]))

    inds 
end

