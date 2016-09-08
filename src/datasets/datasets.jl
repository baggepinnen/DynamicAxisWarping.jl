module Datasets

import BinDeps: unpack_cmd

# UCR datasets #
export ucr_traindata, ucr_testdata, download_ucr

# synthetic datasets #
export fakedata

# place to store the data
const DATAPATH = Pkg.Dir.path()*"/TimeWarp/data"

include("ucr_datasets.jl")
include("fake_datasets.jl")

end # DATASETS module
