

# main dev routine for Mopt.jl
cd(Pkg.dir("MOpt"))

# load the source of the package
include("src/MOpt.jl")

# test a full example
mm = MOpt.serialNormal()

# test a full example
mm = MOpt.BGP_example()


# to develop with tests: run this
include("test/test_MProb.jl")
include("test/test_chains.jl")
include("test/test_BGPchain.jl")
include("test/test_algoBGP.jl")
include("test/test_slices.jl")
include("test/test_sobol.jl")
include("test/test_eval.jl")
include("test/test_objfunc.jl")
include("test/runtests.jl")
		

# to develop in main: run this
# runnign examples:
# MOpt.save(MA,"algo.h5")
	

pb    = ["a" => [1.9,-2,2] , "b" => [-0.9,-1,1] ] 
moms = MOpt.DataFrame(name=["alpha","beta"],value=[0.0,0.0],weight=rand(2))

# run a dummy objective fuction
mprob = @> MOpt.MProb() MOpt.addSampledParam!(pb) MOpt.addMoment!(moms) MOpt.addEvalFunc!(MOpt.Testobj2)
MA = MOpt.MAlgoBGP(mprob,opts)
MOpt.runMOpt!(MA)


# run a bivariate normal objective fuction
moms = MOpt.DataFrame(name=["mu1","mu2"],value=[0.0,0.0],weight=[0.1,0.1])
mprob = @> MOpt.MProb() MOpt.addSampledParam!(pb) MOpt.addMoment!(moms) MOpt.addEvalFunc!(MOpt.objfunc_norm)
MA = MOpt.MAlgoBGP(mprob,opts)
MOpt.runMOpt!(MA)
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
 m = MProb()
 addParam!(m,"c",-1)
 addSampledParam!(m,"a")
 addSampledParam!(m,"b",-0.9,-1,1)
end


