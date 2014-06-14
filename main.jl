

# main dev routine for Mopt.jl

home = ENV["HOME"]
cd("$home/git/MOpt.jl")

# to develop with tests: run this
include("test/test_MProb.jl")
include("test/test_chains.jl")
include("test/test_BGPchain.jl")
include("test/test_algoBGP.jl")
	
# to develop in main: run this
include("src/mopt.jl")



p    = ["a" => 3.1 , "b" => 4.9]
pb   = [ "a" => [0,1] , "b" => [0,1] ]
moms = [
	"alpha" => [ 0.8 , 0.02 ],
	"beta"  => [ 0.8 , 0.02 ],
	"gamma" => [ 0.8 , 0.02 ]
]

mprob = Mopt.MProb(p,pb,Mopt.Testobj,moms)
opts =["N"=>5,"shock_var"=>1.0,"mode"=>"serial","maxiter"=>100,"path"=>".","maxtemp"=>100,"min_shock_sd"=>1.0,"max_shock_sd"=>15.0,"past_iterations"=>30,"min_jumptol"=>0.1,"max_jumptol"=>1.0] 
MA = Mopt.MAlgoBGP(mprob,opts)
Mopt.runMopt(MA)


# workflow
# =========

# step 1: define a MProb
# ----------------------

# get a parameter vector
p = ["a" => 3.1 , "b" => 4.9]
# define params to use with bounds
pb= [ "a" => [0,1.0] , "b" => [0.0,1] ]	

# get some moments
# first entry is moment estimate, second is standard deviation
moms = [
	"alpha" => [ 0.8 , 0.02 ],
	"beta"  => [ 0.8 , 0.02 ],
	"gamma" => [ 0.8 , 0.02 ]
]

# a subset of moments to match
submoms = ["alpha", "beta"]

# call objective
x = Mopt.Testobj(p,moms,submoms)

# Define an Moment Optimization Problem
mprob = Mopt.MProb(p,pb,Mopt.Testobj,moms;moments_subset=submoms)

# show
mprob



# step 2: choose an algorithm
# ----------------------
opts =["N"=>3,"shock_var"=>1.0,"mode"=>"serial","maxiter"=>100,"path"=>"."] 

MA = Mopt.MAlgoRandom(mprob,opts)


# an algo has N chains:
MA.MChains
# step 3: run estimation
# ----------------------
runMopt(algo)




# start iteration
updateChain!(chains,mprob,mprob.initial_value)

for i in 2:chainLength
	p = computeNewGuess(algo,chains,mprob)
	updateChain!(chains,mprob,mprob.initial_value)
end





