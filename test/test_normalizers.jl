using DynamicAxisWarping: advance!, lastlength
using Statistics

n = 10
x = randn(100)
z = ZNormalizer(x,n)

@test length(z) == n
@test size(z) == (n,)
@test ndims(z) == 1
@test lastlength(z) == length(z)

inds = 1:n
@test z == x[inds]
@test mean(z) ≈ mean(x[inds])
@test std(z) ≈ std(x[inds], corrected=false)
@test z.i == 1

@test z[!,1] == (x[1]-mean(z))/std(z)
@test z[!,n] == (x[n]-mean(z))/std(z)

advance!(z)
inds = inds .+ 1

@test z.i == 2
@test z == x[inds]
@test mean(z) ≈ mean(x[inds])
@test std(z) ≈ std(x[inds], corrected=false)

@test z[!,1] == (x[2]-mean(z))/std(z)
@test z[!,n] == (x[n+1]-mean(z))/std(z)



@test_throws BoundsError z[n+1]

for i = 1:91
    advance!(z)
end
@test_throws BoundsError z[n]
