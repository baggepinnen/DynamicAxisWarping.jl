using DynamicAxisWarping, Plots

fs = 100
t  = range(0,stop=1,step=1/fs)
y0 = sin.(2pi .*t)
y1 = sin.(2.5pi .*t)
y  = [y0;y1[2:end]]
q  = [y0;y0[2:end]]

f1 = plot([q y])
f2 = dtwplot(q,y,lc=:green, lw=1)
f3 = matchplot(q,y,ds=3,separation=1)
plot(f1,f2,f3, legend=false, layout=3)
