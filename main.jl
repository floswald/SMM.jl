

# main dev routine for Mopt.jl
cd(Pkg.dir("MOpt"))

# load the source of the package
include("src/MOpt.jl")

# test a full example
mm = MOpt.serialNormal()

# to develop with tests: run this
include("test/test_MProb.jl")
include("test/test_eval.jl")
include("test/test_BGPchain.jl")
include("test/test_algoBGP.jl")
include("test/test_slices.jl")
include("test/test_sobol.jl")
include("test/test_objfunc.jl")
include("test/runtests.jl")
		


