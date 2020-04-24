using DynamicAxisWarping
using DynamicAxisWarping.Datasets
using Plots
gr()

data, labels = ucr_traindata("50words")

plot(
    [fastdtw(data[:, i], data[:, i+10], radius)[1] for radius = 1:100, i = 1:10],
    show = true,
)
xlabel!("radius")
