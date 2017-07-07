using TimeWarp
 #using TimeWarp.WarpPlots
using TimeWarp.Datasets
using PyPlot # use pyplot
plt = PyPlot

# UCI data repository must be downloaded first, run TimeWarp.Datasets.download_ucr()
data, labels = ucr_traindata("MedicalImages");

# c = [ cl==1 ? :blue : :red for cl in class ]'
# plot(y, linecolor=c, legend=false)
nclust = 3
init_centers = TimeWarp.dbaclust_initial_centers(data, nclust, FastDTW(10))
centers, clustids, result = dbaclust(data, nclust, FastDTW(10); init_centers=deepcopy(init_centers), iterations=1)

figure()
plot(data, color="0.75")
plot(centers, color="red")
plot(init_centers, color="green")

if is_linux()
  plt.show()
end
