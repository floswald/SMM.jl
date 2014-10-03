

# main dev routine for Mopt.jl

home = ENV["HOME"]
cd("$home/git/MOpt.jl")

# load the source of the package
include("src/MOpt.jl")

# test a full example
include("src/cluster/examples/example-serial.jl")

# to develop with tests: run this
include("test/test_MProb.jl")
include("test/test_chains.jl")
include("test/test_BGPchain.jl")
include("test/test_algoBGP.jl")
	
include("test/runtests.jl")
		

# to develop in main: run this

# runnign examples:
	
# MOpt.save(MA,"algo.h5")


p    = ["a" => 1.9 , "b" => -0.9]
pb   = [ "a" => [-2,2] , "b" => [-1,1] ]
moms = MOpt.DataFrame(moment=["alpha","beta"],data_value=[0.0,0.0],data_sd=rand(2))
# moms = [
# 	"alpha" => [ 0.8 , 0.02 ],
# 	"beta"  => [ 0.8 , 0.02 ],
# 	"gamma" => [ 0.8 , 0.02 ]
# ]

mprob = MOpt.MProb(p,pb,MOpt.objfunc_norm2,moms)
MOpt.MAlgoBGP(mprob,opts)

x = MOpt.slices(mprob,30)
MOpt.plotSlices(mprob,x[1],x[2])
	
# usign Gadfly
# MOpt.plotSlices(mprob,x,joinpath(pwd(),"slices.pdf"))

## Testing new interface
using MOpt

mprob = @> begin
 MProb()
 addParam!("c",-1)
 addSampledParam!("a",1.9,-2,2)
 addSampledParam!("b",-0.9,-1,1)
end


