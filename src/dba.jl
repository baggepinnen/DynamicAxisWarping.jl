"""
    DBAResult(cost,converged,iterations)

Holds basic results of a DTW Barycenter Averaging (DBA) fit.
"""
type DBAResult
    cost::Float64
    converged::Bool
    iterations::Int
    cost_trace::Vector{Float64}
end

"""
    avgseq, results = dba(sequences, [dist=SqEuclidean()]; kwargs...)

Perfoms DTW Barycenter Averaging (DBA) given a collection of `sequences`
and the current estimate of the average sequence

Example usage:

    x = [1,2,2,3,3,4]
    y = [1,3,4]
    z = [1,2,2,4]
    avg,result = dba([x,y,z])
"""
function dba{T<:AbstractVecOrMat}(
        sequences::AbstractVector{T},
        dist::SemiMetric = SqEuclidean();
        n::Int = 0,
        iterations::Int = 1000,
        rtol::Float64 = 1e-3,
        store_trace::Bool = false
    )

    # initialize dbavg as signal with length n
    if n <= 0
        n = mean([ length(s) in sequences ])
        dbavg = zeros(n)
    end

    # variables storing optimization progress
    converged = false
    iter = 0
    cost = Inf
    cost_trace = Float64[]

    ## main loop ##
    while !converged && iter < iterations
        # do an iteration of dba
        newavg, newcost = dba_iteration(dbavg,sequences,dist)
        iter += 1

        # store history of cost while optimizing (optional)
        store_trace && push!(cost_trace, newcost)

        # check convergence
        if (cost-newcost)/newcost > rtol
            converged = true
        else
            cost = newcost
        end

        # update estimate
        dbavg = newavg
    end

    return dbavg, DBAResult(cost,converged,iter,trace)
end


"""
    newavg, cost = dba_iteration(dbavg, sequences, dist)

Performs one iteration of DTW Barycenter Averaging (DBA) given a collection of
`sequences` and the current estimate of the average sequence, `dbavg`. Returns
an updated estimate, and the cost/loss of the previous estimate
"""
function dba_iteration{T<:AbstractVecOrMat}(
        dbavg::T,
        sequences::AbstractVector{T},
        dist::SemiMetric
    )

    count = zeros(Int, length(dbavg))
    newavg = zeros(Float64, length(dbavg))
    total_cost = 0.0
    
    for seq in sequences
        # time warp signal versus average
        cst, match1, match2 = dtw(dbavg, seq)
        total_cost += cost
        
        # store stats for barycentric average
        for j=1:length(match2)
            count[i1[j]] += 1
            newavg[i1[j]] += seq[i2[j]]
        end
    end

    # compute average and return total cost
    newavg = newavg ./ count
    return newavg, total_cost
end
