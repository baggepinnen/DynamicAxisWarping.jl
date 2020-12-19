export dtwplot


"""
    dtwplot(seq1, seq2, [dist=SqEuclidean()]; transportcost=1, diagonal=false)
    dtwplot(seq1, seq2, D, i1, i2; transportcost=1, diagonal=false)

Given two sequences, perform dynamic time warping and plot
the results. If alignment has already been computed, pass
the indices `i1` and `i2` to make the plot.

`diagonal = true` plots a diagonal marker as visual aid.
"""
dtwplot

function handleargs(seq1, seq2, dist::SemiMetric = SqEuclidean(); kwargs...)
    D = dtw_cost_matrix(seq1, seq2, dist; kwargs...)
    cost, i1, i2 = DynamicAxisWarping.trackback(D)
    seq1, seq2, D, i1, i2
end

function handleargs(seq1, seq2, dist, i2min, i2max; kwargs...)
    D = dtw_cost_matrix(seq1, seq2, dist, i2min, i2max; kwargs...)
    cost, i1, i2 = DynamicAxisWarping.trackback(D)
    seq1, seq2, D, i1, i2
end

handleargs(h; kwargs...) = handleargs(h.args...; kwargs...)

@userplot DTWPlot

@recipe function f(h::DTWPlot; transportcost=1, diagonal=false)
    seq1, seq2, D, i1, i2 = handleargs(h; transportcost=transportcost)

    n1, n2 = lastlength(seq1), lastlength(seq2)

    all = ndims(seq1) ∈ (1,2)
    # set up the subplots
    legend --> false
    link := :both
    grid --> false
    if all
        layout --> @layout [
            a b{0.8w,0.8h}
            _ c
        ]
    else
        layout --> 1
    end

    left_margin --> 0mm
    bottom_margin --> 0mm
    top_margin --> 0mm
    right_margin --> 0mm

    # heatmap
    @series begin
        clims --> (0, 3 * D[end, end])
        seriestype := :heatmap
        formatter --> (z) -> ""
        subplot := (all ? 2 : 1)
        D
    end

    # the rest of the plots are paths

    # main plot
    s1 = @series begin
        seriestype := :path
        linecolor --> :auto
        linewidth --> 3
        subplot := (all ? 2 : 1)
        formatter --> (z) -> ""
        i1, i2
    end

    if diagonal
        m2 = max(n1, n2)
        m1 = min(n1, n2)
        d = m2-m1
        seriestype := :path
        subplot := (all ? 2 : 1)
        imi, ima = radiuslimits(d,seq1, seq2)
        if d == 0
            @series 1:n1
        else
            @series begin
                [imi ima]
            end
        end
    end

    if all
        if ndims(seq1) == 1
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
        else
            # left line plot
            @series begin
                seriestype := :heatmap
                subplot := 1
                seq2'
            end

            # bottom line plot
            @series begin
                seriestype := :heatmap
                subplot := 3
                seq1
            end

        end
    end
end



@userplot MatchPlot

znorm(x) = (x = x.- mean(x); x ./= std(x))
using Statistics
@recipe function f(h::MatchPlot; transportcost=1, separation=2, ds=1)
    x, y, D, i1, i2 = handleargs(h; transportcost=transportcost)
    x,y = znorm.((x,y))
    s1 = x .- separation
    s2 = y .+ separation

    @series begin
        s1
    end
    @series begin
        s2
    end
    @series begin
        primary := false
        linecolor --> :black
        seriesalpha --> 0.2
        i = fill(Inf, 1, length(i1))
        vec([i1'; i2'; i][:,1:ds:end]), vec([s1[i1]'; s2[i2]'; i][:,1:ds:end])
    end
end

@recipe function plot(r::DTWSearchResult)
    title --> "DTW-NN Search result"
    yguide --> "Distance"
    label --> round(r.cost, sigdigits=4)
    if length(r.dists) == 1
        @series begin
            seriestype := :scatter
            makersize --> 15
            markershape --> :x
            group := 1
            r.dists
        end
        @series begin
            seriestype := :hline
            linestyle := :dash
            primary := false
            group := 1
            r.dists
        end
    else
        @series r.dists
    end

    @series begin
        seriestype := :vline
        linestyle := :dash
        linecolor := :black
        primary := false
        [r.loc]
    end

end
