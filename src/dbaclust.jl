"""
    DBAclustResult(centers,clustids,result)

Holds results of a DBAclust run.

"""
mutable struct DBAclustResult{T}
    centers::T
    clustids::Array{Int}
    converged::Bool
    iterations::Int
    dbaresult::DBAResult
end

"""

    optional_threaded(cond, ex)

Execute ex with multi-threading if cond is true,
otherwise execute with one core
"""
macro optional_threaded(cond, ex)
    quote
        if $(esc(cond))
            $(esc(:(Threads.@threads $ex)))
        else
            $(esc(ex))
        end
    end
end


"""
    dbaclust(
        sequences,
        nclust::Int,
        dtwdist::DTWDistance;
        n_init::Int           = 1,
        iterations::Int       = 100,
        inner_iterations::Int = 10,
        rtol::Float64         = 1e-4,
        rtol_inner::Float64   = rtol,
        threaded::Bool        = false,
        show_progress::Bool   = true,
        store_trace::Bool     = true,
        i2min::AbstractVector = [],
        i2max::AbstractVector = [],
    )


# Arguments:
- `nclust`: Number of clsuters
- `n_init`: Number of initialization tries
- `inner_iterations`: Number of iterations in the inner alg.
- `i2min`: Bounds on the warping path
- `threaded`: Use multi-threading 
"""
function dbaclust(
    sequences,
    nclust::Int,
    dtwdist::DTWDistance;
    n_init::Int           = 1,
    iterations::Int       = 100,
    inner_iterations::Int = 10,
    rtol::Float64         = 1e-4,
    rtol_inner::Float64   = rtol,
    threaded::Bool        = false,
    show_progress::Bool   = true,
    store_trace::Bool     = true,
    i2min::AbstractVector = [],
    i2max::AbstractVector = [],
)

    n_init < 1 && throw(ArgumentError("n_init must be greater than zero"))

    T = typeof(sequences)

    results = Array{DBAclustResult{T}}(undef, n_init)
    show_progress && (p = Progress(n_init))

    @optional_threaded threaded for i = 1:n_init
        results[i] = dbaclust_single(
            sequences,
            nclust,
            dtwdist;
            iterations = iterations,
            inner_iterations = inner_iterations,
            rtol = rtol,
            rtol_inner = rtol_inner,
            threaded = threaded,
            show_progress = false,
            store_trace = store_trace,
            i2min = i2min,
            i2max = i2max,
        )
        show_progress && next!(p)
    end
    best = results[1]

    for i = 2:n_init
        if results[i].dbaresult.cost < best.dbaresult.cost
            best = results[i]
        end
    end

    return best
end


"""
    avgseq, results = dbaclust_single(sequences, dist; kwargs...)

Perfoms a single DTW Barycenter Averaging (DBA) given a collection of `sequences`
and the current estimate of the average sequence.

Example usage:

    x = [1,2,2,3,3,4]
    y = [1,3,4]
    z = [1,2,2,4]
    avg,result = dba([x,y,z], DTW(3))
"""
function dbaclust_single(
    sequences::AbstractVector,
    nclust::Int,
    dtwdist::DTWDistance;
    threaded::Bool        = false,
    init_centers::AbstractVector = dbaclust_initial_centers(
        sequences,
        nclust,
        dtwdist;
        threaded
    ),
    iterations::Int       = 100,
    inner_iterations::Int = 10,
    rtol::Float64         = 1e-4,
    rtol_inner::Float64   = rtol,
    show_progress::Bool   = true,
    store_trace::Bool     = true,
    i2min::AbstractVector = [],
    i2max::AbstractVector = [],
)

    T = floattype(eltype(sequences))
    # rename for convienence
    avgs   = init_centers
    N = length(avgs[1])

    # check initial centers have the same length
    if !all(length(a) == N for a in avgs)
        throw(ArgumentError("all initial centers should be the same length"))
    end

    # dimensions
    nseq      = length(sequences)
    maxseqlen = maximum([length(s) for s in sequences])


    # TODO switch to ntuples?
    counts    = [zeros(Int, N) for _ = 1:nclust]
    sums      = [Array{T}(undef,N) for _ = 1:nclust]

    # cluster assignments for each sequence
    clus_asgn = Array{Int}(undef,nseq)
    c         = 0

    # arrays storing path through dtw cost matrix
    i1, i2    = Int[], Int[]

    # variables storing optimization progress
    converged       = false
    iter            = 0
    inner_iter      = 0
    converged_inner = false
    last_cost       = Inf
    total_cost      = 0.0
    cost_trace      = Float64[]
    costs           = Array{Float64}(undef,nseq)

    # main loop ##
    if show_progress
        prog = ProgressMeter.ProgressThresh(rtol, 2)
        ProgressMeter.update!(
            prog,
            Inf;
            showvalues = [
                (:iteration, iter),
                (Symbol("max iteration"), iterations),
                (:cost, total_cost),
            ],
        )
    end#showprogress

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
            if threaded
                # using multi-core
                cluster_dists_      = Array{Float64}(undef, nclust)
                cluster_i1s_        = Array{Vector{Int64}}(undef, nclust)
                cluster_i2s_        = Array{Vector{Int64}}(undef, nclust)

                Threads.@threads for c_ = 1:nclust
                    # if one of the two is empty, use unconstrained window. If both are nonempty, but not the same lenght, distpath will throw error
                    if isempty(i2min) && isempty(i2max)
                        cost, i1_, i2_ = distpath(dtwdist, avgs[c_], seq)
                    else
                        cost, i1_, i2_ = distpath(dtwdist, avgs[c_], seq, i2min, i2max)
                    end
                    cluster_dists_[c_]  = cost
                    cluster_i1s_[c_]    = i1_
                    cluster_i2s_[c_]    = i2_
                end
                cost, c = findmin(cluster_dists_)
                i1      = cluster_i1s_[c]
                i2      = cluster_i2s_[c]
                costs[s] = cost
            else
                # using single-core
                for c_ = 1:nclust
                    # if one of the two is empty, use unconstrained window. If both are nonempty, but not the same lenght, distpath will throw error
                    if isempty(i2min) && isempty(i2max)
                        cost, i1_, i2_ = distpath(dtwdist, avgs[c_], seq)
                    else
                        cost, i1_, i2_ = distpath(dtwdist, avgs[c_], seq, i2min, i2max)
                    end
                    if cost < costs[s]
                        # store cluster, and match indices
                        c           = c_
                        i1          = i1_
                        i2          = i2_
                        costs[s]    = cost
                    end
                end
            end

            # s was assigned to cluster c
            clus_asgn[s] = c
            cnt          = counts[c]
            sm           = sums[c]
            avg          = avgs[c]

            # update stats for barycentric average for
            # the assigned cluster
            for t = 1:length(i2)
                cnt[i1[t]] += 1
                sm[i1[t]]  += seq[i2[t]]
            end
        end

        # if any centers are unused, and reassign them to the sequences
        # with the highest cost
        unused = setdiff(1:nclust, unique(clus_asgn))
        if !isempty(unused)
            # reinitialize centers
            @optional_threaded threaded for c in unused
                avgs[c] = deepcopy(sequences[argmax(costs)])
                @optional_threaded threaded for s = 1:nseq
                    seq = sequences[s]
                    if isempty(i2min) && isempty(i2max)
                        cost, = distpath(dtwdist, avgs[c], seq)
                    else
                        cost, = distpath(dtwdist, avgs[c], seq, i2min, i2max)
                    end
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
        Δ = (last_cost - total_cost) / total_cost
        if Δ < rtol
            converged = true
        else
            last_cost = total_cost
        end

        # update barycenter estimates
        for (a, s, c) in zip(avgs, sums, counts)
            for t = 1:N
                c[t] == 0 && continue
                a[t] = s[t] / c[t]
            end
            # zero out sums and counts for next iteration
            s .= 0
            c .= 0
        end

        # add additional inner dba iterations

        for i = 1:nclust
            seqs            = view(sequences, clus_asgn .== i)
            inner_iter      = 0
            converged_inner = false
            oldcost         = 1.0e100
            while !converged_inner && inner_iter < inner_iterations
                newcost = dba_iteration!(
                    sums[i],
                    avgs[i],
                    counts[i],
                    seqs,
                    dtwdist;
                    i2min = i2min,
                    i2max = i2max,
                )
                copy!(avgs[i], sums[i])
                inner_iter += 1
                δ = (oldcost - newcost) / oldcost
                if δ < rtol_inner
                    converged_inner = true
                else
                    oldcost = newcost
                end
            end
        end

        # update progress bar
        iter += 1
        show_progress && ProgressMeter.update!(
            prog,
            Δ;
            showvalues = [
                (:iteration, iter),
                (Symbol("max iteration"), iterations),
                (:cost, total_cost),
            ],
        )
    end

    return DBAclustResult(
        avgs,
        clus_asgn,
        converged,
        iter,
        DBAResult(total_cost, converged_inner, inner_iter, cost_trace),
    )
end


"""
   dbaclust_initial_centers(sequences, nclust, dtwdist::DTWDistance; threaded::Bool = false)

Uses kmeans++ (but with dtw distance) to initialize the centers
for dba clustering.
"""
function dbaclust_initial_centers(
    sequences::AbstractVector,
    nclust::Int,
    dtwdist::DTWDistance;
    threaded::Bool        = false
)
    # number of sequences in dataset
    nseq          = length(sequences)
    # distance of each datapoint to each center
    dists         = zeros(nclust, nseq)
    # distances to closest center
    min_dists     = zeros(1, nseq)
    # choose a center uniformly at random
    center_ids    = zeros(Int, nclust)
    center_ids[1] = rand(1:nseq)

    # assign the rest of the centers
    for c = 1:(nclust-1)

        # first, compute distances for the previous center
        cent = sequences[center_ids[c]]
        @optional_threaded threaded for i = 1:nseq
            # this distance will be zero
            i == center_ids[c] && continue
            # else, compute dtw distance
            seq = sequences[i]
            dists[c, i], = distpath(dtwdist, seq, cent)
        end

        # for each sequence, find distance to closest center
        minimum!(min_dists, dists[1:c, :])

        min_dists .= abs2.(min_dists)

        # sample the next center
        center_ids[c+1] = sample(1:nseq, Weights(view(min_dists, :)))
    end

    # return list of cluster centers
    return [deepcopy(sequences[c]) for c in center_ids]
end
