
module TestMProb

using FactCheck, DataFrames

# eventually:
# put push!(LOAD_PATH, "/Path/To/Mopt/module/") into your  ~\.juliarc.jl
# using Mopt

# for now:
include("../src/mopt.jl")


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

	pb   = (ASCIIString=>Array{Float64,1})[ "a" => [0,1] , "b" => [0,1] ]
	moms = [
		"alpha" => [ 0.8 , 0.02 ],
		"beta"  => [ 0.8 , 0.02 ],
		"gamma" => [ 0.8 , 0.02 ]
	]
	mprob = Mopt.MProb(p,pb,Mopt.Testobj,moms);
	@fact collect(Mopt.ps_names(mprob)) == collect(keys(pb)) => true
	@fact collect(Mopt.ms_names(mprob)) == collect(keys(moms)) => true

end



end # module 






