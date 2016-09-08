using TimeWarp
using Plots
gr()

# UCI data repository must be downloaded first, run TimeWarp.download_data()
data, labels = TimeWarp.traindata("MedicalImages");

# c = [ cl==1 ? :blue : :red for cl in class ]'
# plot(y, linecolor=c, legend=false)
nclust = 3
init_centers = TimeWarp.dbaclust_initial_centers(data, nclust, FastDTW(10))
centers, clustids, result = dbaclust(data, nclust, FastDTW(10); init_centers=deepcopy(init_centers), iterations=1)

plot(data, line=(:grey), legend=false)
plot!(centers, line=(:red), legend=false)
plot!(init_centers, line=(:green), legend=false)

