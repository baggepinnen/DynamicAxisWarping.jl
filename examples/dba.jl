using TimeWarp
 #using TimeWarp.WarpPlots
using TimeWarp.Datasets

# UCI data repository must be downloaded first, run TimeWarp.Datasets.download_ucr()
data, labels = ucr_traindata("Gun_Point");

# c = [ cl==1 ? :blue : :red for cl in class ]'
# plot(y, linecolor=c, legend=false)

avg, result = dba(data, FastDTW(1))

plot(data)
plot(avg)
