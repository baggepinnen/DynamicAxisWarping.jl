using DynamicAxisWarping, SlidingDistancesBase
using DelimitedFiles, ReverseDiff, DiffResults, Optim
cd(@__DIR__)

data = readdlm("../data/CBF_TRAIN.txt")
c1,c2,c3 = ntuple(i->data[data[:,1] .== i, 2:end], 3)

c1v = copy.(collect(eachrow(c1)))

plot(c1', layout=3, sp=1, legend=false)
plot!(c2', sp=2)
plot!(c3', sp=3)

function dtwmeancost(data, γ)
    function (x)
        mean(data) do d
            soft_dtw_cost(x,d,γ=γ)/lastlength(d)
        end
    end
end

input = mean(c1v) # Our initial guess will be the Euclidean mean
costfun = dtwmeancost(c1v, 1)
costfun(input)

cfg = ReverseDiff.GradientConfig(input)
tape = ReverseDiff.GradientTape(costfun, input)
ctape = ReverseDiff.compile(tape)
results = DiffResults.GradientResult(similar(input))


function fg!(F,G,x)
    if G != nothing
        ReverseDiff.gradient!(results, ctape, x)
        G .= results.derivs[1]
        return results.value
    end
    if F != nothing
        return costfun(x)
    end
end

##
res = Optim.optimize(
    Optim.only_fg!(fg!),
    input,
    BFGS(),
    Optim.Options(
        store_trace       = true,
        show_trace        = true,
        show_every        = 1,
        iterations        = 100,
        allow_f_increases = false,
        time_limit        = 100,
        x_tol             = 1e-3,
        f_tol             = 1e-4,
        g_tol             = 1e-4,
    ),
)

##
using Plots, Plots.PlotMeasures
f1 = plot(c1', lab="", axis=false, legend=:bottom)
plot!(input, l=(4, :red), lab="Euclidean mean")
plot!(res.minimizer, l=(4, :green), lab="Soft-DTW mean")
f2 = plot(c1[1:3,:]', layout=(1,3), legend=false, axis=false, margin=-5mm)
plot(f1,f2, layout=(2,1))
