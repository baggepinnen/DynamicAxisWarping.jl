module Datasets

# synthetic datasets #
export fakedata

# place to store the data
const DATAPATH = joinpath(@__DIR__, "../../data") |> normpath

include("fake_datasets.jl")

end # DATASETS module
