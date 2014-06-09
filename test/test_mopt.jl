
module TestMopt

using FactCheck, DataFrames

# eventually:
# put push!(LOAD_PATH, "/Path/To/Mopt/module/") into your  ~\.juliarc.jl
# using Mopt

# for now:
include("../src/mopt.jl")


# create test data:
# 1) a test objective function
# 2) parameter Dict
# 3) parameters to estimate Dict
# 4) moments dict
# 5) moments subset Array{String}



# Choose Algorithm and configure it

# malgo = MAlgoRandom(mprob,3)
# malgo["maxiter"]   = 100
# malgo["save_freq"] = 5
# malgo["mode"]      = "serial"


# TESTING MProb
# ==============

# test default constructor type
p    = ["a" => 3.1 , "b" => 4.9]
pb   = [ "a" => [0,1] , "b" => [0,1] ]
moms = [
	"alpha" => [ 0.8 , 0.02 ],
	"beta"  => [ 0.8 , 0.02 ],
	"gamma" => [ 0.8 , 0.02 ]
]

facts("Testing the MProb constructor") do

	context("default constructor") do

		mprob = Mopt.MProb(p,pb,Mopt.Testobj,moms)
		@fact isa(mprob, Mopt.MProb) => true

	end

	context("constructor throws errors") do

		# get an inexisting moment to subset
		@fact_throws Mopt.MProb(p,pb,Mopt.Testobj,moms; moments_subset=["alpha","epsilon"]);

		# get some wrong moments
		moms = [
			"alpha" => [ 0.8 , 0.02 ],
			"beta"  => [ 0.8 , 0.02 ],
			"gamma" => [ 0.8  ]
		]
		@fact_throws Mopt.MProb(p,pb,Mopt.Testobj,moms);

		# i'm giving a parameter "c" that is not in initial_value
		pb= [ "a" => [0,1] , "c" => [0,1] ]
		moms = [
			"alpha" => [ 0.8 , 0.02 ],
			"beta"  => [ 0.8 , 0.02 ],
			"gamma" => [ 0.8 , 0.02 ]
		]
		@fact_throws Mopt.MProb(p,pb,Mopt.Testobj,moms);

		# get some wrong bounds
		pb= [ "a" => [1,0] , "b" => [0,1] ]
		@fact_throws Mopt.MProb(p,pb,Mopt.Testobj,moms);
	end

end


facts("testing MProb methods") do

	pb   = [ "a" => [0,1] , "b" => [0,1] ]
	moms = [
		"alpha" => [ 0.8 , 0.02 ],
		"beta"  => [ 0.8 , 0.02 ],
		"gamma" => [ 0.8 , 0.02 ]
	]
	mprob = Mopt.MProb(p,pb,Mopt.Testobj,moms);
	@fact collect(Mopt.ps_names(mprob)) == collect(keys(pb)) => true
	@fact collect(Mopt.ms_names(mprob)) == collect(keys(moms)) => true

end



# TESTING Chains
# ==============


p    = ["a" => 3.1 , "b" => 4.9]
pb   = [ "a" => [0,1] , "b" => [0,1] ]
moms = [
	"alpha" => [ 0.8 , 0.02 ],
	"beta"  => [ 0.8 , 0.02 ],
	"gamma" => [ 0.8 , 0.02 ]
]

facts("Testing Chains constructor") do
	
	mprob = Mopt.MProb(p,pb,Mopt.Testobj,moms)
	L = 9
	chain = Mopt.Chain(mprob,L)

	@fact chain.i => 0 

	context("length of members") do

		# test that all member except i are L long
		@fact length(chain.evals) => L
		@fact length(chain.accept) => L
		for nm in Mopt.ps_names(mprob)
			@fact length(chain.parameters[nm]) => L
		end
		for nm in Mopt.ms_names(mprob)
			@fact length(chain.moments[nm]) => L
		end
	end

	context("names of param and moments dicts") do

		@fact collect(keys(chain.parameters)) == collect(keys(mprob.initial_value)) => true
		@fact collect(keys(chain.moments)) == collect(keys(mprob.moments)) => true

	end

end

facts("testing Chains methods") do
	
	context("test appendEval!") do

		mprob = Mopt.MProb(p,pb,Mopt.Testobj,moms)
		v = Mopt.Testobj(p,moms,["alpha","beta","gamma"])
		L = 9
		chain = Mopt.Chain(mprob,L)

		# verify values are zero:
		@fact all(chain.evals .== 0.0) => true
		@fact chain.i => 0 

		# set i to 1 to test this:
		chain.i = 1

		# update chain with v
		Mopt.appendEval!(chain,v)

		# verify new values on chain
		@fact chain.i => 1 
		@fact chain.evals[1] => v["value"]
		for nm in Mopt.ps_names(mprob)
			@fact chain.parameters[nm][1] => v["params"][nm]
		end
		for nm in Mopt.ms_names(mprob)
			@fact chain.moments[nm][1] => v["moments"][nm]
		end

	end



end





end # module TestMopt

# test dummy functions





