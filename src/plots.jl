export dtwplot


"""
    dtwplot(seq1, seq2, [dist=SqEuclidean()])
    dtwplot(seq1, seq2, D, i1, i2)

Given two sequences, perform dynamic time warping and plot
the results. If alignment has already been computed, pass
the indices `i1` and `i2` to make the plot.
"""
function handleargs(seq1, seq2, dist::SemiMetric = SqEuclidean(); kwargs...)
    D = dtw_cost_matrix(seq1, seq2, dist; kwargs...)
    cost, i1, i2 = DynamicAxisWarping.trackback(D)
    seq1, seq2, D, i1, i2
end

handleargs(h; kwargs...) = handleargs(h.args...; kwargs...)

@userplot DTWPlot

@recipe function f(h::DTWPlot; transportcost=1)
    seq1, seq2, D, i1, i2 = handleargs(h; transportcost=transportcost)

    n1, n2 = lastlength(seq1), lastlength(seq2)

    all = ndims(seq1) == 1
    # set up the subplots
    legend --> false
    link := :both
    grid --> false
    if all
        layout --> @layout [
            a b{0.9w,0.9h}
            _ c
        ]
    else
        layout --> 1
    end

    left_margin --> 0mm
    bottom_margin --> 0mm
    top_margin --> 0mm
    right_margin --> 0mm
    clims --> (0, 3 * D[end, end])

    # heatmap
    @series begin
        seriestype := :heatmap
        formatter --> (z) -> ""
        subplot := (all ? 2 : 1)
        D
    end

    # the rest of the plots are paths
    seriestype := :path
    linecolor --> RGB(0, 0, 0)

    # main plot
    s1 = @series begin
        linewidth --> 3
        subplot := (all ? 2 : 1)
        formatter --> (z) -> ""
        i1, i2
    end

    if all
        # left line plot
        @series begin
            subplot := 1
            seq2, 1:n2
        end

        # bottom line plot
        @series begin
            subplot := 3
            1:n1, seq1
        end
    end
end
