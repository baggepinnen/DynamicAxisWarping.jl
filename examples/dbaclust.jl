using TimeWarp
using Plots
pyplot()

# UCI data repository must be downloaded first, run TimeWarp.download_data()
data, labels = TimeWarp.traindata("gun_point");

# c = [ cl==1 ? :blue : :red for cl in class ]'
# plot(y, linecolor=c, legend=false)

centers, clustids, result = dbaclust(data, 3)

plot(data, linecolor=:grey, legend=false)
for avg in centers
    plot!(avg, linecolor=:red, legend=false, line=(2))
end

