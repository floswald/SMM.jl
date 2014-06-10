


module TestAlgo

using FactCheck, DataFrames

include("../src/mopt.jl")


p    = ["a" => 3.1 , "b" => 4.9]
pb   = [ "a" => [0,1] , "b" => [0,1] ]
moms = [
	"alpha" => [ 0.8 , 0.02 ],
	"beta"  => [ 0.8 , 0.02 ],
	"gamma" => [ 0.8 , 0.02 ]
]

mprob = Mopt.MProb(p,pb,Mopt.Testobj,moms)
L = 9
chain = Mopt.Chain(mprob,L)

facts("testing MAlgoRandom Constructor") do

	context("checking members") do

		opts =["N"=>3,"shock_var"=>1.0,"mode"=>"serial","maxiter"=>100,"path"=>"."] 

		MA = Mopt.MAlgoRandom(mprob,opts)

		@fact isa(MA,Mopt.MAlgo) => true
		@fact isa(MA,Mopt.MAlgoRandom) => true
		@fact isa(MA.m,Mopt.MProb) => true

		@fact MA.i => 0
		@fact MA.chains.n => opts["N"]

	end

	context("checking getters/setters on MAlgo") do

		opts =["N"=>3,"shock_var"=>1.0,"mode"=>"serial","maxiter"=>100,"path"=>"."] 
		MA = Mopt.MAlgoRandom(mprob,opts)

		# getters
		@fact MA["N"] => opts["N"]
		@fact MA["mode"] => opts["mode"]

		# setter
		MA["newOption"] = "hi"
		@fact MA["newOption"] => "hi"

	end




end







end # module

# context("test updateChain!()") do

	# 	mprob = Mopt.MProb(p,pb,Mopt.Testobj,moms)
	# 	v = Mopt.Testobj(p,moms,["alpha","beta","gamma"])
	# 	L = 9
	# 	chain = Mopt.Chain(mprob,L)

	# 	# verify values are zero:
	# 	@fact all(chain.evals .== 0.0) => true
	# 	@fact chain.i => 0 

	# 	# call updateChain!
	# 	Mopt.updateChain!(chain,mprob,p)

	# 	# verify new values on chain
	# 	@fact chain.evals[1] => v["value"]
	# 	for nm in Mopt.ps_names(mprob)
	# 		@fact chain.parameters[nm][1] => v["params"][nm]
	# 	end
	# 	for nm in Mopt.ms_names(mprob)
	# 		@fact chain.moments[nm][1] => v["moments"][nm]
	# 	end

	# end
