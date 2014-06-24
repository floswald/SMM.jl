

# main dev routine for Mopt.jl

home = ENV["HOME"]
cd("$home/git/MOpt.jl")

# to develop with tests: run this
include("test/runtests.jl")

include("test/test_MProb.jl")
include("test/test_chains.jl")
include("test/test_BGPchain.jl")
include("test/test_algoBGP.jl")
	
	
# to develop in main: run this
include("src/MOpt.jl")

# runnign examples:
banan=false
include("src/examples.jl")
	
MOpt.save(MA,"algo.h5")


p    = ["a" => 3.1 , "b" => 4.9]
pb   = [ "a" => [0,1] , "b" => [0,1] ]
moms = DataFrame(moment=["alpha","beta","gamma"],data_value=[0.8,0.7,0.5],data_sd=rand(3))
# moms = [
# 	"alpha" => [ 0.8 , 0.02 ],
# 	"beta"  => [ 0.8 , 0.02 ],
# 	"gamma" => [ 0.8 , 0.02 ]
# ]

MOpt.Testobj(p,moms,["alpha","beta","gamma"])
mprob = MOpt.MProb(p,pb,MOpt.Testobj,moms)
# opts =["N"=>20,"shock_var"=>1.0,"mode"=>"serial","maxiter"=>100,"path"=>".","maxtemp"=>100,"min_shock_sd"=>1.0,"max_shock_sd"=>1.0,"past_iterations"=>30,"min_jumptol"=>0.1,"max_jumptol"=>1.0] 
# MA = Mopt.MAlgoBGP(mprob,opts)
# Mopt.runMopt(MA)





