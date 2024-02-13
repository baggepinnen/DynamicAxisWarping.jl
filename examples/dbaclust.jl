using DynamicAxisWarping
#using DynamicAxisWarping.WarpPlots
using DynamicAxisWarping.Datasets

# UCI data repository must be downloaded first, then extract and prepare the data set labeled "MedicalImages"
# data, labels = ... # MedicalImages

# c = [ cl==1 ? :blue : :red for cl in class ]'
# plot(y, linecolor=c, legend=false)
nclust = 3
init_centers = DynamicAxisWarping.dbaclust_initial_centers(data, nclust, FastDTW(10))
centers, clustids, result =
  dbaclust(data, nclust, FastDTW(10); init_centers = deepcopy(init_centers), iterations = 1)

plot(data)
plot(centers)
plot(init_centers)
