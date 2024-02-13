using DynamicAxisWarping
using DynamicAxisWarping.Datasets
using Plots
gr()

# UCI data repository must be downloaded first, then extract and prepare the data set labeled "50words"
# data, labels = ... # 50words

# data, labels = # 50words data...

plot(
    [fastdtw(data[:, i], data[:, i+10], radius)[1] for radius = 1:100, i = 1:10],
    show = true,
)
xlabel!("radius")
