using DynamicAxisWarping
 #using DynamicAxisWarping.WarpPlots
using DynamicAxisWarping.Datasets

# UCI data repository must be downloaded first, then extract and prepare the data set labeled "Gun_Point"
# data, labels = ... # Gun_Point

# c = [ cl==1 ? :blue : :red for cl in class ]'
# plot(y, linecolor=c, legend=false)

avg, result = dba(data, FastDTW(1))

plot(data)
plot(avg)
