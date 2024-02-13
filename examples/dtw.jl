using DynamicAxisWarping
using DynamicAxisWarping.Datasets

# UCI data repository must be downloaded first, then extract and prepare the data set labeled "Gun_Point"
# data = ... # Gun_Point

y1 = data[1][:,1]
y2 = data[1][:,2]

#dtwplot(y1,y2,clim=[0,1])
dtwplot(y1,y2)
savefig("dtw.png")
