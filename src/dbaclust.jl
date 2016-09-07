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
function dbaclust{N,T}(
        sequences::AbstractVector{Sequence{N,T}},
        nclust::Int,
        dist::SemiMetric = SqEuclidean();
        centers::AbstractVector{Sequence{N,T}} = dbaclust_initial_centers(sequences, nclust, dist),
        dbalen::Int = 0,
        iterations::Int = 100,
        inner_iterations::Int = 0,
        rtol::Float64 = 1e-4,
        store_trace::Bool = true
    )
    
    # rename for convienences
    avgs = centers

    # dimensions
    nseq = length(sequences)

    # initialize dbavg as signal with length dbalen
    if dbalen <= 0
        dbalen = round(Int,mean([ length(s) for s in sequences ]))
    end

    maxseqlen = maximum([ length(s) for s in sequences ])
    
    # TODO switch to ntuples?
    counts = [ zeros(Int,dbalen) for _ in 1:nclust ]
    sums = [ Sequence(Array(T,dbalen)) for _ in 1:nclust ]

    # cluster assignments for each sequence
    clus_asgn = Array(Int, nseq)
    c = 0

    # arrays storing path through dtw cost matrix
    i1,i2 = Int[],Int[]

    # variables storing optimization progress
    converged = false
    iter = 0
    last_cost = Inf
    total_cost = 0.0
    cost_trace = Float64[]
    costs = Array(Float64, nseq)

    # main loop ##
    prog = ProgressMeter.ProgressThresh(rtol)
    while !converged && iter < iterations
        # first, update cluster assignments based on nearest
        # centroid (measured by dtw distance). Keep track of
        # total cost (sum of all distances to centers).
        total_cost = 0.0
        for s = 1:nseq
            # process sequence s
            seq = sequences[s]

            # find cluster assignment for s
            costs[s] = Inf
            for c_ = 1:nclust
                cost, i1_, i2_ = dtw(avgs[c_], seq, dist)
                if cost < costs[s]
                    # store cluster, and match indices
                    c = c_
                    i1 = i1_
                    i2 = i2_
                    costs[s] = cost
                end
            end

            # s was assigned to cluster c
            clus_asgn[s] = c
            cnt = counts[c]
            sm = sums[c]
            avg = avgs[c]

            # update stats for barycentric average for
            # the assigned cluster
            for t=1:length(i2)
                cnt[i1[t]] += 1
                sm[i1[t]] += seq[i2[t]]
            end
        end

        # if any centers are unused, and reassign them to the sequences
        # with the highest cost
        unused = setdiff(1:nclust, unique(clus_asgn))
        if ~isempty(unused)
            # reinitialize centers
            for c in unused
                avgs[c] = deepcopy(sequences[indmax(costs)])
                for s = 1:nseq
                    cost, = dtw(avgs[c], seq, dist)
                    if costs[s] > cost
                        costs[s] = cost
                    end
                end
            end
            # we need to reassign clusters, start iteration over
            continue
        end

        # store history of cost while optimizing (optional)
        total_cost = sum(costs)
        store_trace && push!(cost_trace, total_cost)

        # check convergence
        Δ = (last_cost-total_cost)/total_cost
        if Δ < rtol
            converged = true
        else
            last_cost = total_cost
        end

        # update barycenter estimates
        for (a,s,c) in zip(avgs,sums,counts)
            for t = 1:dbalen
                c[t] == 0 && continue
                a[t] = s[t]/c[t]
            end
            # zero out sums and counts for next iteration
            scale!(s,0)
            scale!(c,0)
        end

        # add additional inner dba iterations
        for (i,a) in enumerate(avgs)
            for inner_iter = 1:inner_iterations
                a = dba_iteration(a,view(sequences, clus_asgn .== i),dist)
            end
        end

        # update progress bar
        iter += 1
        ProgressMeter.update!(prog, Δ; showvalues =[(:iteration,iter),
                                                 (Symbol("max iteration"),iterations),
                                                 (:cost,total_cost)])
    end

    return avgs, clus_asgn, DBAResult(total_cost,converged,iter,cost_trace)
end

# Wrapper for AbstractArray of one-dimensional time series.
function dbaclust( s::AbstractArray, args...; kwargs... )
    dbaclust(_sequentize(s), args...; kwargs...)
end

"""
   dbaclust_initial_centers(sequences, nclust, dist)

Uses kmeans++ (but with dtw distance) to initialize the centers
for dba clustering.
"""
function dbaclust_initial_centers{N,T}(
        sequences::AbstractVector{Sequence{N,T}},
        nclust::Int,
        dist::SemiMetric = SqEuclidean()
    )

    # number of sequences in dataset
    nseq = length(sequences)
    
    # distance of each datapoint to each center
    dists = zeros(nclust, nseq)

    # distances to closest center
    min_dists = zeros(1, nseq)

    # choose a center uniformly at random
    center_ids = zeros(Int,nclust)
    center_ids[1] = rand(1:nseq)

    # assign the rest of the centers
    for c = 1:(nclust-1)

        # first, compute distances for the previous center
        cent = sequences[center_ids[c]]
        for (i,seq) in enumerate(sequences)
            # this distance will be zero
            i == center_ids[c] && continue
            # else, compute dtw distance
            dists[c,i], = dtw(seq, cent, dist)
        end

        # for each sequence, find distance to closest center
        minimum!(min_dists, dists[1:c,:])

        # square distances
        map!((x)->x^2, min_dists)

        # sample the next center
        center_ids[c+1] = sample(1:nseq, WeightVec(view(min_dists,:)))
    end

    # return list of cluster centers
    return [ deepcopy(sequences[c]) for c in center_ids ]
end

# Wrapper for AbstractArray of one-dimensional time series.
function dbaclust_initial_centers(
        s::AbstractArray,
        nclust::Int,
        dist::SemiMetric = SqEuclidean()
    )
    dbaclust_initial_centers(_sequentize(s), nclust, dist)
end