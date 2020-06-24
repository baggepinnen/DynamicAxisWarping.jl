using DynamicAxisWarping, Plots

ts = range(0, stop=4π, length=128)
x = LinearInterpolation(sin.(ts) .+ 0.1 .* randn.())
y = LinearInterpolation(sin.(1.1 .* ts))

t = range(0, stop = 1, length=100)
orig_plot = plot(t, x, label="x(t)")
plot!(t, y, label="y(t)")
xlabel!("Time t")
ylabel!("Amplitude")
title!("Original signals")


cost, ϕ, ψ = gdtw(x, y)
gdtw_plot = plot(t, x ∘ ϕ, label="x∘ϕ(t)")
plot!(t, y ∘ ψ, label="y∘ψ(t)")
xlabel!("Time t")
ylabel!("Amplitude")
title!("GDTW")

dtw_cost, t_dtw_x, t_dtw_y = dtw(x.x, y.x)
dtw_plot = plot(eachindex(t_dtw_x) ./ length(t_dtw_x), x.x[t_dtw_x], label="x(t_x)")
plot!(eachindex(t_dtw_x) ./ length(t_dtw_x), y.x[t_dtw_y], label="y(t_y)")
xlabel!("Time")
ylabel!("Amplitude")
title!("DTW")


warp_plot = plot(t, ϕ, label = "ϕ(t)")
plot!(t, ψ, label = "ψ(t)")
rescale_t_x = (t_dtw_x .- first(t_dtw_x)) ./ last(t_dtw_x)
rescale_t_y = (t_dtw_y .- first(t_dtw_y)) ./ last(t_dtw_y)
xlabel!("Time")
ylabel!("Rescaled time")
title!("Comparison of warpings")
plot!(t, LinearInterpolation(rescale_t_x), label="t_x")
plot!(t, LinearInterpolation(rescale_t_y), label="t_y")


plot(gdtw_plot, dtw_plot, orig_plot, warp_plot, layout = (2, 2), legend = :outertopright)

