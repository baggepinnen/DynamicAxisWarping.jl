using TimeWarp
using TimeWarp.WarpPlots
gr() # use pyplot

# UCI data repository must be downloaded first, run TimeWarp.download_data()
data, labels = TimeWarp.traindata("gun_point");

# c = [ cl==1 ? :blue : :red for cl in class ]'
# plot(y, linecolor=c, legend=false)

avg, result = dba(data, FastDTW(1))

plot(data, linecolor=:grey, legend=false)
plot!(avg, linecolor=:red, legend=false, line=(2))
