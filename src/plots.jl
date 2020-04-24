export dtwplot


"""
    dtwplot(seq1, seq2, [dist=SqEuclidean()])
    dtwplot(seq1, seq2, D, i1, i2)

Given two sequences, perform dynamic time warping and plot
the results. If alignment has already been computed, pass
the indices `i1` and `i2` to make the plot.
"""
function dtwplot(
        seq1,
        seq2,
        dist::SemiMetric=SqEuclidean();
        kwargs...
    )
    D = dtw_cost_matrix(seq1, seq2, dist; kwargs...)
    cost,i1,i2 = DynamicAxisWarp.trackback(D)
    dtwplot(seq1, seq2, D, i1, i2)
end

@userplot DTWPlot

@recipe function f(h::DTWPlot)
    seq1, seq2, D, i1, i2 = h.args
    n1, n2 = length(seq1), length(seq2)

    # set up the subplots
    legend --> false
    link := :both
    grid --> false
    layout --> @layout [
        a  b{0.9w,0.9h}
        _  c
    ]

    left_margin --> 0mm
    bottom_margin --> 0mm
    top_margin --> 0mm
    right_margin --> 0mm
    clim --> (0,3*D[end,end])

    # heatmap
    @series begin
        seriestype := :heatmap
        formatter --> (z)->""
        subplot := 2
        D
    end

    # the rest of the plots are paths
    seriestype := :path
    linecolor --> RGB(0,0,0)

    # main plot
    @series begin
        linewidth --> 3
        subplot := 2
        formatter --> (z)->""
        i1, i2
    end

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
