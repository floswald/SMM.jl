
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

	# get model moments as a dict
	mdict = {i => mm[i,:model] for i in keys(mom) }

	# subset to required moments only
	mm = mm[findin(mm[:name],whichmom),:]

	# compute distance
	v = sum((mm[:data] - mm[:model]).^2)

	# status
	status = 1

	# time out
	t0 = time() - t0

	# return a dict
	ret = ["value" => v, "param" => x, "time" => t0, "status" => status, "moments" => mdict]
	return ret

end

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
		moms = [
			"alpha" => [ 0.8 , 0.02 ],
			"beta"  => [ 0.8 , 0.02 ],
			"gamma" => [ 0.8 , 0.02 ]
		]
		@fact_throws Mopt.MProb(p,pb,Testobj,moms);

		# get some wrong bounds
		pb= [ "a" => [1,0] , "b" => [0,1] ]
		@fact_throws Mopt.MProb(p,pb,Testobj,moms);
	end

end


facts("testing MProb methods") do

	pb   = [ "a" => [0,1] , "b" => [0,1] ]
	moms = [
		"alpha" => [ 0.8 , 0.02 ],
		"beta"  => [ 0.8 , 0.02 ],
		"gamma" => [ 0.8 , 0.02 ]
	]
	mprob = Mopt.MProb(p,pb,Testobj,moms);
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
	
	mprob = Mopt.MProb(p,pb,Testobj,moms)
	L = 9
	chain = Mopt.Chain(mprob,L)

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





end # module TestMopt

# test dummy functions





