
module TestMProb

using FactCheck, MOpt




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
moms = DataFrame(moment=["alpha","beta","gamma"],data_value=[0.8,0.7,0.5],data_sd=rand(3))
# moms = [
# 	"alpha" => [ 0.8 , 0.02 ],
# 	"beta"  => [ 0.8 , 0.02 ],
# 	"gamma" => [ 0.8 , 0.02 ]
# ]

facts("Testing the MProb constructor") do

	context("default constructor") do

		mprob = MProb(p,pb,Testobj,moms)
		@fact isa(mprob, MProb) => true

	end

	context("constructor throws errors") do

		# get an inexisting moment to subset
		@fact_throws MProb(p,pb,Testobj,moms; moments_subset=["alpha","epsilon"]);

		# get some wrong moments
		# moms = [
		# 	"alpha" => [ 0.8 , 0.02 ],
		# 	"beta"  => [ 0.8 , 0.02 ],
		# 	"gamma" => [ 0.8  ]
		# ]
		# @fact_throws MProb(p,pb,Testobj,moms);

		# i'm giving a parameter "c" that is not in initial_value
		pb= [ "a" => [0,1] , "c" => [0,1] ]
		# moms = [
		# 	"alpha" => [ 0.8 , 0.02 ],
		# 	"beta"  => [ 0.8 , 0.02 ],
		# 	"gamma" => [ 0.8 , 0.02 ]
		# ]
		@fact_throws MProb(p,pb,Testobj,moms);

		# get some wrong bounds
		pb= [ "a" => [1,0] , "b" => [0,1] ]
		@fact_throws MProb(p,pb,Testobj,moms);
	end

end


facts("testing MProb methods") do

	pb   = (ASCIIString=>Array{Float64,1})[ "a" => [0,1] , "b" => [0,1] ]
	# moms = [
	# 	"alpha" => [ 0.8 , 0.02 ],
	# 	"beta"  => [ 0.8 , 0.02 ],
	# 	"gamma" => [ 0.8 , 0.02 ]
	# ]
	mprob = MProb(p,pb,Testobj,moms);
	@fact collect(MOpt.ps_names(mprob)) == collect(keys(pb)) => true
	@fact MOpt.ms_names(mprob) == mprob.moments[:moment] => true

	@fact eltype(mprob.p2sample_sym)==Symbol => true

end


end # module 






