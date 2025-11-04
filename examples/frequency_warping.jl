using DynamicAxisWarping, DSP, LPVSpectral, Distances, SpectralDistances

# Create two signals
N = 48_000
g(x,N) = exp(-10*(x-N/2)^2/N^2)
t = 1:N
f = range(0.01, stop=1, length=N)
y = sin.(t .* f) .* g.(t, N)
y1 = [y; 0y] .+ 0.1 .* randn.()
y2 = [0y; y] .+ 0.1 .* randn.()

plot([y1 y2])
#
M1,M2 = melspectrogram.((y1,y2), 2048)
plot(plot(M1), plot(M2))
# Calculate the Dynamic Time-Warping cost between them (I do this in √ domain since that appears to produce the nicest looking plots)
@btime dtw(.√(M1.power), .√(M2.power))[1] evals=1 samples=5
# Calculating this distance is rather fast.
# Visualize the coupling matrix. The green line should cross diagonally (since the two signals are the same) starting halfway along the axes (since one signal is shifted by 50% of the length copared to the other).
dtwplot(.√(M1.power), .√(M2.power), linecolor=:green)
#

# The code above used the squared Euclidean distance between spectra for each time point. We can replace that distance with a transport-based distance instead.
# We thus combine dynamic time warping with a transport-based cost along the frequency axis
n,m = size(M1.power)
dist = DiscreteGridTransportDistance(Cityblock(), eltype(M1.power), n, n)
dtw(.√(M1.power), .√(M2.power), dist)[1]
#
dtwplot(.√(M1.power), .√(M2.power), dist, linecolor=:green)
# We now see that the valley through which the path can move is slightly wider. This happens since there is now an additional degree of freedom due to transport along the frequency axis. Since there is no need for that transport in this case, the resulting path remains almost the same.
# We can create the need for some transport along the frequency axis

f1 = range(0.01, stop=1, length=N)
f2 = range(0.1, stop=1, length=N).^2
y10 = sin.(t .* f1) .* g.(t, N)
y20 = sin.(t .* f2) .* g.(t, N)
y1 = [y10; 0y10] .+ 0.1 .* randn.()
y2 = [0y20; y20] .+ 0.1 .* randn.()
M1,M2 = melspectrogram.((y1,y2), 2048)
plot(plot(M1), plot(M2))
# We now redo the same stuff as above
dtw(.√(M1.power), .√(M2.power))[1]
#
dtwplot(.√(M1.power), .√(M2.power), linecolor=:green)
#
dtw(.√(M1.power), .√(M2.power), dist)[1]
#
dtwplot(.√(M1.power), .√(M2.power), dist, linecolor=:green)
# We expect to see less time warping (straighter green line) now when it's allowed to fudge the distance a bit by shuffling mass along the frequency axis.
# How long time does it take to compute such a distance? (on a laptop, expect about 3x improvement on desktop)
@btime dtw($(M1.power), $(M2.power), $dist)
@btime dtw($(M1.power), $(M2.power), radiuslimits(50,m,m)..., $dist)
# The size of the spectrogram
size(M1.power)
# The complexity is almost linear in the number of frequency bins, but quadratic in the number of time steps.
