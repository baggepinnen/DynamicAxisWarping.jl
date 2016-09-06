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
function dba{N,T}(
        sequences::AbstractVector{Sequence{N,T}},
        nclust::Int = 1,
        dist::SemiMetric = SqEuclidean();
        dbalen::Int = 0,
        iterations::Int = 1000,
        rtol::Float64 = 1e-5,
        store_trace::Bool = false
    )

    # dimensions
    nseq = length(sequences)

    # initialize dbavg as signal with length dbalen
    if dbalen <= 0
        dbalen = round(Int,mean([ length(s) for s in sequences ]))
    end

    maxseqlen = maximum([ length(s) for s in sequences ])
    
    avgs = ntuple((_)->Sequence(randn(dbalen)/sqrt(dbalen)), nclust)
    counts = ntuple((_)->zeros(Int,dbalen), nclust) 
    sums = ntuple((_)->Sequence(Array(T,dbalen), nclust)

    # cluster assignments for each sequence
    clus_asgn = Array(Int, nseq)

    # arrays storing path through dtw cost matrix
    i1,i2 = Int[],Int[]

    # variables storing optimization progress
    converged = false
    iter = 0
    last_cost = Inf
    total_cost = 0.0
    cost_trace = Float64[]

    # main loop ##
    prog = ProgressMeter.ProgressThresh(rtol)
    while !converged && iter < iterations
        
        # first, update cluster assignments based on nearest
        # centroid (measured by dtw distance). Keep track of
        # total cost (sum of all distances to centers).
        for s = 1:nseq
            # process sequence s
            seq = sequences[s]

            # find cluster assignment for s
            min_cost = Inf
            for c_ = 1:nclust
                cost, i1_, i2_ = dtw(avgs[c_], seq, dist)
                if cost < min_cost
                    # store cluster, and match indices
                    c = c_
                    i1 = i1_
                    i2 = i2_
                    min_cost = cost
                end
            end

            # s was assigned to cluster c
            clus_asgn[s] = c
            cnt = counts[c]
            sm = sums[c]
            avg = avgs[c]

            # the contribution of s to the total cost
            total_cost += min_cost

            # update stats for barycentric average for
            # the assigned cluster
            for t=1:length(i2)
                cnt[i1[t]] += 1
                sm[i1[t]] += seq[i2[t]]
                # TODO: double check direction here
            end
        end

        # store history of cost while optimizing (optional)
        store_trace && push!(cost_trace, newcost)

        # check convergence
        Δ = (last_cost-total_cost)/total_cost
        if Δ < rtol
            converged = true
        else
            last_cost = total_cost
        end

        # update barycenter estimates
        for a,s,c in zip(avgs,sums,counts)
            for t = 1:dbalen
                a[t] = s[t]/c[t]
            end
        end

        # update progress bar
        iter += 1
        ProgressMeter.update!(prog, Δ; showvalues =[(:iteration,iter),
                                                 (Symbol("max iteration"),iterations),
                                                 (:cost,cost)])
    end

    return avgs, clus_asgn, DBAResult(cost,converged,iter,cost_trace)
end
