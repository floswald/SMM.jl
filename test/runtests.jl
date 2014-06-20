
module test_MOpt

using FactCheck
using MOpt

include("test_chains.jl")
include("test_MProb.jl")
include("test_algoBGP.jl")
include("test_BGPchain.jl")

exitstatus()

end


