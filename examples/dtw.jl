using TimeWarp
using TimeWarp.Datasets

# UCI data repository must be downloaded first, run TimeWarp.download_data()
data = ucr_traindata("Gun_Point");

y1 = data[1][:,1]
y2 = data[1][:,2]

#dtwplot(y1,y2,clim=[0,1])
dtwplot(y1,y2)
savefig("dtw.png")
