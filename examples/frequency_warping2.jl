# # Finding patterns in signals using Optimal Transport and Dynamic Time Warping

# In this example we will compare different distances when used as the inner metric in a dynamic time-warping of spectrograms.
# Setup
using DynamicAxisWarping, SpectralDistances, LPVSpectral, DSP, Plots, Statistics, BenchmarkTools, Distances, Random, ThreadTools, AlphaStableDistributions
theme(:default)
plotly();
Random.seed!(0);

# We start by creating a short patter, the *query*, we then create a long time series that contains a quite similar sound, but where the frequencies of the chirp are slightly higher. This is a realistic scenario when doing, e.g., acoustic detection. Different individuals within the same animal species might have similar calls, but slightly different pitch etc. We also add some alpha sub-gaussian noise in order to make the problem a bit harder.

N       = 48_000 ÷ 2
nfft    = 4096
f1      = range(0.01, stop = 1, length = N)
f2      = range(1, stop = 0.1, length = N) .^ 2
g(x, N) = exp(-10 * (x - N / 2)^2 / N^2)
t       = 1:N
y10     = (sin.(t .* f1) .+ sin.(t .* (f2 .+ 0.4))) .* g.(t, N)
y20     = (sin.(t .* (f1 .+ 0.5)) .+ sin.(t .* (f2 .+ 0.5))) .* g.(t, N)
q       = Float32.(y10 .+ 0.2 .* randn.())
y       = Float32.([0y20; 0y20; y20; 0y20; 0y20] .+ 1 .* randn.())
y     .+= v1(rand(AlphaSubGaussian(n=length(y))));
# y = repeat(y, 40) #src

Q = melspectrogram(q, nfft, window = hanning, fmin = 400, fs = 48_000)
Y = melspectrogram(y, nfft, window = hanning, fmin = 400, fs = 48_000)
plot(plot(Q, title = "Query"), plot(Y, title = "Data"), link = :both, layout = (2, 1))
# We now see if we can detect the pattern using DTW with the standard squared Euclidean distance
rad  = 10 # This is the maximum allowed warping radius
dist = SqEuclidean()
w    = DTWWorkspace(sqrt.(Q.power), dist, rad)
res  = dtwnn(w, sqrt.(Y.power))
plot(Y);
vline!([Y.time[res.loc]], l = (4, :blue), primary=false)
# That did probably not go well at all! The line indicates where the smallest distance to the pattern was, i.e., where the nearest neighbor search thinks the *onset* of the pattern is.
#---
# Let's do the same with a transport-based distance
n, m = size(Q.power)
dist = DiscreteGridTransportDistance(Cityblock(), n, n)
w    = DTWWorkspace(sqrt.(Q.power), dist, rad)
res  = dtwnn(w, sqrt.(Y.power))
plot(Y);
vline!([Y.time[res.loc]], l = (4, :blue), primary=false)
# The line should now be a better indication of the onset of the pattern
#---
# Next, we plot the cost function for each time shift to see how they behave
using DynamicAxisWarping: lastlength
function naive(a, b, dist = SqEuclidean(), r = 7)
    dists = map(1:lastlength(b)-lastlength(a)) do i
        dtw_cost(a, b[!, i:i+lastlength(a)-1], dist, r)
    end
end

plot()
for dist in [DiscreteGridTransportDistance(Cityblock(), n, n), SqEuclidean()]
    @time res = naive(Q.power, Y.power, dist, rad)
    res .-= minimum(res)
    res ./= maximum(res)
    plot!(res, lab = string(typeof(dist).name))
end
fig1 = vline!([2 * size(Y.power, 2) / 5], l = (4, :blue), primary=false)
# The Euclidean distance is rather terrible for this task. The transport based distance seems to be doing the right thing.
#---
# unfortunately, calculating this using the transport-based distance takes quite a while
@btime dtwnn($w, $(sqrt.(Y.power)), prune_endpoints=false)
# while that might not seem like much, the data we're searching through only corresponds to this many seconds, if the sample rate is 48kHz
length(y) / 48000
#---
# We are thus intersted in distances that are faster to compute. We'll explore the OTRD distance, which performs optimal transport between the roots of a linear-system representation of the signal.
fm  = TimeWindow(inner = LS(na = 8, λ=1e-3), n = nfft, noverlap = nfft÷2)
Qm  = fm(q) |> change_precision(Float64)
Ym  = fm(y) |> change_precision(Float64);
# Let's see how the world looks like through the lens of autoregressive models
plot(Qm, LPVSpectral.mel_to_hz(Q.mels)./48000, rad=false)
#---
dist  = OptimalTransportRootDistance(p = 1, β = 0.5)
distf = (x, y) -> evaluate(dist, x, y, tol = 1e-3)
res   = naive(Qm.models, Ym.models, distf, rad)
res .-= minimum(res)
res ./= maximum(res)
plot!(fig1, res, lab = "OTRD")
# This should do an okay job, but is it fast?
#---
w = DTWWorkspace(Qm.models, distf, rad)
@btime dtwnn($w, $(Ym.models), prune_endpoints=false)
# About the same, but it could be much faster if the allocations were taken care of.
# Let's speed up even further, using the root distance (RD)
Qe     = Float32.(reduce(hcat, embedding.(Qm.models))) # This extracts the roots of the linear system into a real vector
Ye     = Float32.(reduce(hcat, embedding.(Ym.models)))
dist   = SqEuclidean();
# dist   = EuclideanRootDistance(p=1, weight=e->sqrt.(residueweight(e))) #src
res    = naive(Qe, Ye, dist, rad)
res  .-= minimum(res)
res  ./= maximum(res)
plot!(res, lab = "RD")
# This should do an okay job, but is it fast?
#---
w = DTWWorkspace(Qe, dist, rad)
dtwnn(w, (Ye), prune_endpoints=false)
@btime dtwnn($w, $(Ye), prune_endpoints=false)
# That should be fast very, but there are still some room for improvements through early-termination heuristics.
GC.gc(true); GC.gc(true) #src

# ## Notes on futher performance improvements
# Searching for the nearest neighbor using DTW tends to become relatively cheaper the longer the data searched through is. This is due to pruning heauristics that abort the search early based on the smallest distance found so far.
# The main thing to watch out for is the length (in time) of the query sequence. Using a larger `nfft` effectively reduces the time resolution which greatly reduces the time-complexity, 𝒪(t²) if pruning heuristics are not effective.
# Let's see how long time it takes to search through a much longer time series (note that estimating the models also takes quite a while, but this is 𝒪(t) expensive and can be easily parallelized so I don't count that for now)
yl  = Float32.(repeat(y, 400) .+ 0.2 .* randn.())
println("Length of singal: ", length(yl)/48_000/60, " minutes, at 48kHz")
Yml = fm(yl) |> change_precision(Float32);
Yel = reduce(hcat, embedding.(Yml.models))
sleep(2) #src
@btime dtwnn($w, $(Yel), prune_endpoints=true)


# ### Future work
# The mass normalization for spectrograms has to be improved. Currently, each spectrum (per time point) is normalized to sum to one. A perhaps better strategy is to normalize the entire spectrogram and then use unbalanced mass transport per time point.

GC.gc(true); GC.gc(true) #src

#---
# # Transport in time
# We can also consider using optimal transport along the time axis instead of DTW. Let's do the same stuff as above for the TimeDistance from SpectralDistances.jl
function naive_time(a, b, dist)
    n = length(a)
    dists = map(1:length(b.models)-length(a.models)) do i
        Y = TimeVaryingAR(b.models[i:i+length(a.models)-1])
        dist(a, Y, tol=1e-3, check_interval=2)
    end
end
dist = TimeDistance(inner=OptimalTransportRootDistance(p=1, β=0.5),tp=1,c=0.1)

res2  = naive_time(Qm, Ym, dist)
res   = SpectralDistances.distance_profile(dist, Qm, Ym)
res .-= minimum(res)
res ./= maximum(res)
plot!(fig1, res, lab = "Time transport")
#---
@btime naive_time($Qm, $Ym, $dist);
@btime SpectralDistances.distance_profile($dist, $(change_precision(Float32,Qm)), $(change_precision(Float32,Ym)), check_interval=5)
# It appears to be surprisingly competitive, both in terms of accuracy and performance, even though it's evaluated using a rather naive method. Since the solver used to find the `TimeDistance` internally solves the *dual* problem, a lower bound of the objective function is always available for free during the optimization. This lower bound could be used to terminate the optimization early if it rises above the smallest distance found so far.

@test mean(abs2, SpectralDistances.distance_profile(dist, Qm, Ym) - naive_time(Qm, Ym, dist)) < 1e-4
