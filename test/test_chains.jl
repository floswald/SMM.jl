module TestMChain

using FactCheck, DataFrames

include("../src/mopt.jl")


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

facts("testing Chain methods") do
	
	context("test appendEval!()") do

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



end





end # module 






