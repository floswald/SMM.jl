
module TestMopt

using FactCheck

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


# test constructor# define a Test objective function
function Testobj(x::Dict,mom::Dict,whichmom::Array{ASCIIString,1})

	t0 = time()

	# perform computations of objective function
	# create DataFrame to substract model from data:
	mm = DataFrame(name=collect(keys(mom)),data=collect(values(mom)),model=[x["a"] + x["b"] + i for i=1:nrow(mm)])

	# subset to required moments only
	mm = mm[findin(mm[:name],whichmom),:]

	# compute distance
	v = sum((mm[:data].-mm[:model]).^2)

	# status
	status = 1

	# time out
	t0 = time() - t0

	# return a dict
	ret = ["value" => v, "param" => x, "time" => t0, "status" => status, "moments" => mm]
	return ret

end

# Choose Algorithm and configure it

# malgo = MAlgoRandom(mprob,3)
# malgo["maxiter"]   = 100
# malgo["save_freq"] = 5
# malgo["mode"]      = "serial"


# TESTING
# =======

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

		mprob = Mopt.MProb(p,pb,Testobj,moms)
		@fact isa(mprob, Mopt.MProb) => true

	end

	context("constructor throws errors") do

		# get an inexisting moment to subset
		@fact_throws Mopt.MProb(p,pb,Testobj,moms; moments_subset=["alpha","epsilon"]);

		# get some wrong moments
		moms = [
			"alpha" => [ 0.8 , 0.02 ],
			"beta"  => [ 0.8 , 0.02 ],
			"gamma" => [ 0.8  ]
		]
		@fact_throws Mopt.MProb(p,pb,Testobj,moms);

		# i'm giving a parameter "c" that is not in initial_value
		pb= [ "a" => [0,1] , "c" => [0,1] ]
		@fact_throws Mopt.MProb(p,pb,Testobj,moms);

		# get some wrong bounds
		pb= [ "a" => [1,0] , "b" => [0,1] ]
		@fact_throws Mopt.MProb(p,pb,Testobj,moms);
	end

end






end # module TestMopt

# test dummy functions





