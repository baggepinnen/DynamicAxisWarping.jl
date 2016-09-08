function fakedata(name::AbstractString, args...)
    
    # gaussians with separation
    name == "gaussian" && return fakedata_gaussian(args...)

    # no dataset found
    error("dataset name not recognized.")
end

function fakedata_gaussian(
        pts_per_clust::Int=10,
        nclust::Int=2,
        xmin=0.,
        xmax=nclust*7.,
        nx=round(Int,(xmax-xmin)*10),
        σ=1.0,
        amin=1.,
        amax=2.
    )
    x = linspace(xmin,xmax,nx)
    μ = linspace(xmin,xmax,nclust+2)[2:end-1]
    a = linspace(xmin,xmax,pts_per_clust+2)[2:end-1]
    
    npts = pts_per_clust*nclust
    data = zeros(nx,npts)
    labels = zeros(Int,npts)

    i = 1
    for c = 1:nclust
        for n = 1:pts_per_clust
            data[:,i] = gauss_func(a[n],μ[c],σ,x)
            labels[i] = c
            i += 1
        end
    end

    return data, labels
end

gauss_func(a,μ,σ,x) = a*exp( (-(x-μ).^2) / (2σ^2) )

