
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
moms = DataFrame(name=["alpha","beta","gamma"],value=[0.8,0.7,0.5],weight=rand(3))
# moms = [
# 	"alpha" => [ 0.8 , 0.02 ],
# 	"beta"  => [ 0.8 , 0.02 ],
# 	"gamma" => [ 0.8 , 0.02 ]
# ]

facts("Testing the MProb constructor") do

	context("default constructor") do

		mprob = MProb()
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

	pb   = [ 
		"a" => [0.1, 0, 1] , 
		"b" => [0.2, 0, 1] ]
	# moms = [
	# 	"alpha" => [ 0.8 , 0.02 ],
	# 	"beta"  => [ 0.8 , 0.02 ],
	# 	"gamma" => [ 0.8 , 0.02 ]
	# ]	
	mprob = MProb();
	addSampledParam!(mprob,pb)	
	@fact sort(collect(MOpt.ps_names(mprob))) == sort({:a,:b}) => true

	mprob = MProb();
	addSampledParam!(mprob,"a",0.1,0,1)
	addSampledParam!(mprob,"b",0.1,0,1)
	addMoment(mprob,moms)


	@fact sort(collect(MOpt.ps_names(mprob))) == sort({:a,:b}) => true
	@fact sort(collect(MOpt.ms_names(mprob))) == sort({:alpha,:beta,:gamma}) => true

end


end # module 






