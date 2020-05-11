
function MatrixProfile.matrix_profile(T, m::Int, dist::DTWDistance; showprogress=true)
    n   = lastlength(T)
    l   = n-m+1
    r   = dist.method.radius
    # n > 2m+1 || throw(ArgumentError("Window length too long, maximum length is $((n+1)÷2)"))
    P    = distance_profile(dist, getwindow(T,m,1), T)
    P[1:r] .= typemax(floattype(P))
    D    = similar(P)
    I    = Vector{Int}(undef, l)
    P[1], I[1] = findmin(P)
    prog = Progress((l - 1) ÷ 5, dt=1, desc="Matrix profile", barglyphs = BarGlyphs("[=> ]"), color=:blue)
    bsf = typemax(floattype(T))
    @inbounds for i = 2:l
        Ti = getwindow(T,m,i)
        res = dtwnn(Ti, T, dist.dist, r, transportcost=dist.method.transportcost, avoid=i-r:i+r)
        I[i] = res.loc
        P[i] = res.cost
        showprogress && i % 5 == 0 && next!(prog)
    end
    MatrixProfile.Profile(T, P, I, m, nothing)
end
