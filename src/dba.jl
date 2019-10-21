"""
    DBAResult(cost,converged,iterations,cost_trace)

Holds results of a DTW Barycenter Averaging (DBA) fit.
"""
mutable struct DBAResult
    cost::Float64
    converged::Bool
    iterations::Int
    cost_trace::Vector{Float64}
end

"""
    avgseq, results = dba(sequences, [dist=SqEuclidean()]; kwargs...)

Perfoms DTW Barycenter Averaging (DBA) given a collection of `sequences`
and the current estimate of the average sequence.

Example usage:

    x = [1,2,2,3,3,4]
    y = [1,3,4]
    z = [1,2,2,4]
    avg,result = dba([x,y,z])
"""
function dba(
        sequences::AbstractVector{T},
        method::DTWMethod,
        dist::SemiMetric = SqEuclidean();
        init_center::T = rand(sequences),
        iterations::Int = 1000,
        rtol::Float64 = 1e-5,
        store_trace::Bool = false,
        show_progress::Bool = true,
        i2min::AbstractVector =[],
        i2max::AbstractVector = []
    ) where {T<:Sequence}

    # method for computing dtw
    dtwdist = DTWDistance(method,dist)

    # initialize dbavg as a random sample from the dataset
    nseq = length(sequences)
    dbavg = deepcopy(init_center)

    # storage for each iteration
    newavg = Sequence(zeros(length(dbavg)))
    counts = zeros(Int, length(dbavg))

    # variables storing optimization progress
    converged = false
    iter = 0
    cost,newcost = Inf,Inf
    cost_trace = Float64[]

    # display optimization progress
    if show_progress
        p = ProgressMeter.ProgressThresh(rtol)
    end

    ## main loop ##
    while !converged && iter < iterations

        # do an iteration of dba
        newcost = dba_iteration!(newavg, dbavg, counts, sequences, dtwdist; i2min=i2min,i2max=i2max)
        iter += 1

        # store history of cost while optimizing (optional)
        store_trace && push!(cost_trace, newcost)

        # check convergence
        Δ = (cost-newcost)/newcost
        if Δ < rtol
            converged = true
        else
            # update estimate
            cost = newcost
            dbavg = deepcopy(newavg)
        end

        # update progress bar
        if show_progress
            ProgressMeter.update!(p, Δ; showvalues =[(:iteration,iter),
                                                 (Symbol("max iteration"),iterations),
                                                 (:cost,cost)])
        end
    end

    return newavg, DBAResult(newcost,converged,iter,cost_trace)
end


"""
    newavg, cost = dba_iteration(dbavg, sequences, dist)

Performs one iteration of DTW Barycenter Averaging (DBA) given a collection of
`sequences` and the current estimate of the average sequence, `dbavg`. Returns
an updated estimate, and the cost/loss of the previous estimate
"""
function dba_iteration!(
        newavg::T,
        oldavg::T,
        counts::Array{Int,1},
        sequences::AbstractVector{T},
        d::DTWDistance;
        i2min::AbstractVector=[],
        i2max::AbstractVector=[]
    ) where {T<:Sequence}

    # sum of dtw dist of all sequences to center
    total_cost = 0.0

    # store stats for barycenter averages
    scale!(counts,0)
    scale!(newavg,0)

    # main ploop
    for seq in sequences
        # time warp signal versus average
        # if one of the two is empty, use unconstrained window. If both are nonempty, but not the same lenght, distpath will throw error
        if isempty(i2min) && isempty(i2max)
          cost, i1, i2 = distpath(d, oldavg, seq)
        else
          cost, i1, i2 = distpath(d, oldavg, seq, i2min, i2max)
        end
        total_cost += cost

        # store stats for barycentric average
        for j=1:length(i2)
            counts[i1[j]] += 1
            newavg[i1[j]] += seq[i2[j]]
        end
    end

    # compute average and return total cost
    for i in eachindex(newavg)
        newavg[i] = newavg[i]/counts[i]
    end

    return total_cost
end

# Wrapper for AbstractArray of one-dimensional time series.
function dba( s::AbstractArray, args...; kwargs... )
    dba(_sequentize(s), args...; kwargs...)
end

# weirdly enought, this works for the dtw_dba_miniexample,  but does not work for dtw_dbaclust
#@generated function _sequentize{T,N}(s::AbstractArray{T,N})
#    :( Sequence[ Sequence(@ncall($N, view, s, n-> n==$N ? i : Colon())) for i = 1:size(s,2) ] )
#end


# This function replaces the compile-time function which doesn't work in v0.6
# Compile-time function fixed, generated functions to not support untyped generators
function _sequentize(s::AbstractArray)
    [Sequence(s[:,i]) for i=1:size(s,2)]
end
