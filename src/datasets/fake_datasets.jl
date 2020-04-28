"""
    fakedata_gaussian(pts_per_clust::Int = 10, nclust::Int = 2, xmin = 0.0, xmax = nclust * 7.0, nx = round(Int, (xmax - xmin) * 10), σ = 1.0, amin = 1.0, amax = 2.0)

DOCSTRING

# Arguments:
- `pts_per_clust`: DESCRIPTION
- `nclust`: DESCRIPTION
- `xmin`: DESCRIPTION
- `xmax`: DESCRIPTION
- `nx`: DESCRIPTION
- `σ`: DESCRIPTION
- `amin`: DESCRIPTION
- `amax`: DESCRIPTION
"""
function fakedata_gaussian(
    pts_per_clust::Int = 10,
    nclust::Int        = 2,
    xmin               = 0.0,
    xmax               = nclust * 7.0,
    nx                 = round(Int, (xmax - xmin) * 10),
    σ                  = 1.0,
    amin               = 1.0,
    amax               = 2.0,
)
    x = LinRange(xmin, xmax, nx)
    μ = LinRange(xmin, xmax, nclust + 2)[2:end-1]
    a = LinRange(xmin, xmax, pts_per_clust + 2)[2:end-1]

    npts = pts_per_clust * nclust
    data = zeros(nx, npts)
    labels = zeros(Int, npts)

    i = 1
    for c = 1:nclust
        for n = 1:pts_per_clust
            data[:, i] = gauss_func.(a[n], μ[c], σ, x)
            labels[i] = c
            i += 1
        end
    end

    return data, labels
end

gauss_func(a, μ, σ, x) = a * exp((-(x - μ) .^ 2) / (2 * σ^2))
