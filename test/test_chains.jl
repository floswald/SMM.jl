module TestChain

using FactCheck,DataFrames
include("../src/MOpt.jl")



# TESTING Chains
# ==============
p    = ["a" => 3.1 , "b" => 4.9]
pb   = [ "a" => [0,1] , "b" => [0,1] ]
moms = MOpt.DataFrame(moment=["alpha","beta","gamma"],data_value=[0.8,0.7,0.5],data_sd=rand(3))



facts("Testing Default Chains constructor") do
	
	mprob = MOpt.MProb(p,pb,MOpt.Testobj,moms)
	L = 9
	chain = MOpt.Chain(mprob,L)

	@fact chain.i => 0 

	context("length of members") do

		# test that all member except i are L long
		@fact length(chain.infos[:evals])  => L
		@fact length(chain.infos[:accept]) => L
		for nm in MOpt.ps_names(mprob)
			@fact length(chain.parameters[symbol(nm)]) => L
		end
		for nm in MOpt.ms_names(mprob)
			@fact length(chain.moments[symbol(nm)]) => L
		end
	end
end

# facts("testing MChain constructor") do
	

# 	mprob = MProb(p,pb,Testobj,moms)
# 	L = 9
# 	chain = Chain(mprob,L)

# 	n = 14
# 	chains = MChain(n,mprob,L)

# 	@fact isa(chains,MChain) => true
# 	@fact isa(chains.chains,Array) => true
# 	@fact length(chains.chains) => n
# 	for ch in 1:n
# 		@fact length(chains.chains[ch].infos["evals"]) => L
# 	end
# end

facts("testing Chain/MChain methods") do

	context("testing getindex(chain,i)") do

		mprob = MOpt.MProb(p,pb,MOpt.Testobj,moms)
		v = MOpt.Testobj(p,moms,["alpha","beta","gamma"])
		L = 9	# length of chain
		chain = MOpt.Chain(mprob,L)

		i = rand(1:L)

		x = MOpt.allstats(chain,i)

		@fact isa(x,MOpt.DataFrame) => true
		@fact nrow(x) => 1


	end

	context("testing getindex(chain,i::Range)") do

		mprob = MOpt.MProb(p,pb,MOpt.Testobj,moms)
		v     = MOpt.Testobj(p,moms,["alpha","beta","gamma"])
		L     = 9	# length of chain
		chain = MOpt.Chain(mprob,L)

		i = 3:L
		x = MOpt.allstats(chain,i)
		
		@fact isa(x,MOpt.DataFrame) => true
		@fact nrow(x) => length(i)

	end
	
	context("test appendEval!(chain)") do

		mprob = MOpt.MProb(p,pb,MOpt.Testobj,moms)
		v = MOpt.Testobj(p,moms,["alpha","beta","gamma"])
		for (k,va) in v["params"]
			v["params"][k] = va + 1.1
		end

		L = 9
		chain = MOpt.Chain(mprob,L)

		# set i to 1 to test this:
		chain.i = 1

		# verify values are zero:
		@fact all(chain.infos[chain.i,:evals] == 0.0) => true
		# verify params are zero
		for ic in 1:ncol(chain.parameters)
			@fact chain.parameters[chain.i,ic][1] => 0.0
		end

		# update chain with v
		MOpt.appendEval!(chain,v["value"],v["params"],v["moments"],true,1,rand())

		# verify new values on chain
		@fact chain.infos[:evals][1] => v["value"]
		@fact chain.infos[:accept][1] => true
		for nm in MOpt.ps_names(mprob)
			@fact chain.parameters[chain.i,symbol(nm)][1] => v["params"][nm]
		end
			@fact chain.moments[chain.i,chain.moments_nms] => v["moments"]
	end

	context("test appendEval!(chain) if par is a dataframe") do

		mprob = MOpt.MProb(p,pb,MOpt.Testobj,moms)
		v = MOpt.Testobj(p,moms,["alpha","beta","gamma"])
		for (k,va) in v["params"]
			v["params"][k] = va + 1.1
		end

		L = 9
		chain = MOpt.Chain(mprob,L)

		# set i to 1 to test this:
		chain.i = 1

		# verify values are zero:
		@fact all(chain.infos[chain.i,:evals] == 0.0) => true
		# verify params are zero
		for ic in 1:ncol(chain.parameters)
			@fact chain.parameters[chain.i,ic][1] => 0.0
			# and set random
		end

		oldpar = MOpt.parameters(chain,chain.i)
		for i in chain.params_nms
			oldpar[i] = rand()
		end

		# update chain with v
		MOpt.appendEval!(chain,v["value"],oldpar,v["moments"],true,1,rand())

		# verify new values on chain
		@fact chain.infos[:evals][1] => v["value"]
		@fact chain.infos[:accept][1] => true
		for nm in MOpt.ps_names(mprob)
			@fact chain.parameters[chain.i,symbol(nm)][1] - oldpar[symbol(nm)][1] < 1e-6 => true
		end
			@fact chain.moments[chain.i,chain.moments_nms] => v["moments"]
	end

end


facts("testing collectFields and fillinFields functions") do
	
	mprob = MOpt.MProb(p,pb,MOpt.Testobj,moms)
	chain = MOpt.Chain(mprob,10)

	context("collectFields") do

		df= MOpt.parameters(chain,1)
		@fact isa(df,DataFrame) => true
		@fact df[:a] => [0.0]
		@fact df[:b] => [0.0]

	end

	

end



end # module 






