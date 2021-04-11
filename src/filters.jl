function imfilter(A, kern)
    size(kern,1) == size(kern,2) || throw(ArgumentError("Only square kernels supported, you provided $(size(kern))."))
    isodd(size(kern,1)) || throw(ArgumentError("Only odd-sized kernels supported, you provided $(size(kern))."))
    A2 = LocalFilters.convolve(A, kern)
    pad = size(kern, 1) ÷ 2
    @views A2[1:pad, :] .= A[1:pad, :]
    @views A2[end-pad+1:end, :] .= A[end-pad+1:end, :]
    @views A2[:, 1:pad] .= A[:, 1:pad]
    @views A2[:, end-pad+1:end] .= A[:, end-pad+1:end]
    A2
end


function gaussian2(n)
    t = range(-2,stop=2, length=n)
    k = exp.(-t.^2)
    K = k*k'
    K ./= sum(K)
end

function gaussian(n)
    t = range(-2,stop=2, length=n)
    K = [exp(-(x-y)^2) for x in t, y in t]
    K ./= sum(K)
end