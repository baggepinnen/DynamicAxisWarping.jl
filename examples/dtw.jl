using TimeWarp
using TimeWarp.WarpPlots
pyplot()

# UCI data repository must be downloaded first, run TimeWarp.download_data()
data = TimeWarp.traindata("gun_point");

y1 = vec(data[1,2:end])
y2 = vec(data[2,2:end])

#dtwplot(y1,y2,clim=[0,1])
dtwplot(y1,y2)
savefig("dtw.png")
