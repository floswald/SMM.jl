

# using MOpt 
module MOptTests

using MOpt, FactCheck

include("test_MProb.jl")
include("test_chains.jl")
include("test_BGPchain.jl")
include("test_algoBGP.jl")
include("test_slices.jl")
# include("test_sobol.jl")
include("test_Eval.jl")
include("test_objfunc.jl")
FactCheck.exitstatus()

end



