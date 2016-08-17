"""
    newavg = dba_iteration(dbavg,sequences)

Performs one iteration of DTW Barycenter Averaging (DBA) given a collection of
`sequences` and the current estimate of the average sequence, `dbavg`. Returns
an updated estimate.
"""
function dba_iteration{T<:AbstractVecOrMat}(
        dbavg::T,
        sequences::AbstractVector{T}
    )

    n = length(sequences)
    count = zeros(Int, length(dbavg))
    sumcoords = zeros(Float64, length(dbavg))
    
    for seq in sequences
        cost, match1, match2 = dtw(dbavg, seq)
        for j=1:length(match2)
            count[i1[j]] += 1
            sumcoords[i1[j]] += seq[i2[j]]
        end
        println("Compared $i to the standard")
    end
    sumcoords = sumcoords ./ count
end
