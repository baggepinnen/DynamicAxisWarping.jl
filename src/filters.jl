function imfilter(A::AbstractMatrix, kern)
    out = similar(A)
    @avx for J in CartesianIndices(out)
        tmp = zero(eltype(out))
        for I ∈ CartesianIndices(kern)
            tmp += A[I + J] * kern[I]
        end
        out[J] = tmp
    end
    out
end
