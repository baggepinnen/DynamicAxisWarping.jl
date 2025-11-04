export dtwplot


"""
    dtwplot(seq1, seq2, [dist=SqEuclidean()]; transportcost=1, diagonal=false)
    dtwplot(seq1, seq2, D, i1, i2; transportcost=1, diagonal=false)
    dtwplot(seq1, seq2, D; transportcost=1, diagonal=false)

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

function handleargs(seq1::AbstractArray, seq2::AbstractArray, D::AbstractMatrix; kwargs...)
    cost, i1, i2 = DynamicAxisWarping.trackback(D)
    seq1, seq2, D, i1, i2
end

function handleargs(seq1::AbstractArray, seq2::AbstractArray, D::AbstractMatrix, i1::AbstractVector, i2::AbstractVector; kwargs...)
    seq1, seq2, D, i1, i2
end

function handleargs(seq1, seq2, dist, i2min, i2max; kwargs...)
    D = dtw_cost_matrix(seq1, seq2, i2min, i2max, dist; kwargs...)
    cost, i1, i2 = DynamicAxisWarping.trackback(D)
    seq1, seq2, D, i1, i2
end

handleargs(h; kwargs...) = handleargs(h.args...; kwargs...)

@userplot DTWPlot

@recipe function f(h::DTWPlot; transportcost=1, diagonal=false, postprocess=nothing)
    seq1, seq2, D, i1, i2 = handleargs(h; transportcost=transportcost, postprocess=postprocess)

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

    plots_id = Base.PkgId(Base.UUID("91a5bcdd-55d7-5caf-9e0b-520d859cae80"), "Plots")
    Plots = Base.loaded_modules[plots_id]

    left_margin --> 0Plots.mm
    bottom_margin --> 0Plots.mm
    top_margin --> 0Plots.mm
    right_margin --> 0Plots.mm

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
@recipe function f(h::MatchPlot; transportcost=1, separation=2, ds=1, postprocess=nothing)
    x, y, D, i1, i2 = handleargs(h; transportcost=transportcost, postprocess=postprocess)
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

"""
    matchplot2

Like `matchplot`, but works with 2d and 3d signals.

### Inputs

#### Positional args
For available positional argument patterns, see `DynamicAxisWarping.handleargs`

Generally, it can accept a `x` and `y` signal to warp and optionally `dtw_cost_matrix` inputs, or `x` and `y` signal plus `dtw` or `dtw_cost_matrix` outputs (skipping the warp step).

#### Keyword args
- `transportcost` -- see `dtw_cost_matrix`
- `separation` -- extra separation/padding added between the signals in in ℜⁿ  
- `showindex` -- Whether to add an axis in the plot for the index/"time" axis (appends to the last dimension)

"""
matchplot2

@userplot MatchPlot2
znorm2(x) = (x = x.- mean(x,dims=2); x ./= std(x,dims=2))
@recipe function f(h::MatchPlot2; transportcost=1, separation=0.5, ds=1,
                   postprocess=nothing, showindex=false, normalize=true)

    x, y, D, i1, i2 = DynamicAxisWarping.handleargs(h;
                                                    transportcost=transportcost,
                                                    postprocess=postprocess)
    x,y = normalize ? znorm2.((x,y)) : (x,y)
    if showindex
        x = [x[:,i1]; i1[:,1]']
        y = [y[:,i2]; i2[:,1]']
    else
        x = x[:,i1]
        y = y[:,i2]
    end
    s1 = x .- separation
    s2 = y .+ separation

    @series begin
        (collect(eachrow(s1))...,)
    end
    @series begin
        (collect(eachrow(s2))...,)
    end
    @series begin
        primary := false
        linecolor --> :black
        seriesalpha --> 0.2
        s1, s2 = s1[:,1:ds:end], s2[:,1:ds:end]
        # Concatonate along 3rd dim
        s3 = cat(s1,s2; dims=3)
        s3 = eachrow.(eachslice(s3, dims=2))
        [(rows_of_slices...,) for rows_of_slices in s3]
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
