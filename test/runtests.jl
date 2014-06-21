

# using MOpt can't make that work!
include("../src/MOpt.jl")
using FactCheck

include("test_chains.jl")
include("test_MProb.jl")
include("test_algoBGP.jl")
include("test_BGPchain.jl")

exitstatus()



