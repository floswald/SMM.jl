module TestChain

using FactCheck, MOpt



# TESTING Chains
# ==============
p    = ["a" => 3.1 , "b" => 4.9]
pb   = [ "a" => [0,1] , "b" => [0,1] ]
moms = [
	"alpha" => [ 0.8 , 0.02 ],
	"beta"  => [ 0.8 , 0.02 ],
	"gamma" => [ 0.8 , 0.02 ]
]



facts("Testing Default Chains constructor") do
	
	mprob = MProb(p,pb,Testobj,moms)
	L = 9
	chain = Chain(mprob,L)

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

		mprob = MProb(p,pb,Testobj,moms)
		v = Testobj(p,moms,["alpha","beta","gamma"])
		L = 9	# length of chain
		chain = Chain(mprob,L)

		i = rand(1:L)

		x = allstats(chain,i)

		@fact isa(x,DataFrame) => true
		@fact nrow(x) => 1


	end

	context("testing getindex(chain,i::Range)") do

		mprob = MProb(p,pb,Testobj,moms)
		v     = Testobj(p,moms,["alpha","beta","gamma"])
		L     = 9	# length of chain
		chain = Chain(mprob,L)

		i = 3:L
		x = allstats(chain,i)
		
		@fact isa(x,DataFrame) => true
		@fact nrow(x) => length(i)

	end
	
	context("test appendEval!(chain)") do

		mprob = MProb(p,pb,Testobj,moms)
		v = Testobj(p,moms,["alpha","beta","gamma"])
		L = 9
		chain = Chain(mprob,L)

		# set i to 1 to test this:
		chain.i = 1

		# verify values are zero:
		@fact all(chain.infos[chain.i,:evals] == 0.0) => true


		# update chain with v
		MOpt.appendEval!(chain,v,true,1,rand())

		# verify new values on chain
		@fact chain.infos[:evals][1] => v["value"]
		@fact chain.infos[:accept][1] => true
		for nm in MOpt.ps_names(mprob)
			@fact chain.parameters[chain.i,symbol(nm)][1] => v["params"][nm]
		end
		for nm in MOpt.ms_names(mprob)
			@fact chain.moments[chain.i,symbol(nm)][1] => v["moments"][nm]
		end
	end

	
	context("testing updateIter(MChain") do

		mprob = MProb(p,pb,Testobj,moms)
		L = 9
		n = 17
		MC = [Chain(mprob,L) for i=1:n]

		for ix = 1:n
			@fact MC[ix].i => 0
		end

		MOpt.updateIterChain!(MC)
		for ix = 1:n
			@fact MC[ix].i => 1
		end

	end


end


facts("testing collectFields and fillinFields functions") do
	
	mprob = MProb(p,pb,Testobj,moms)
	chain = Chain(mprob,10)

	context("collectFields") do

		df= parameters(chain,1)
		@fact isa(df,DataFrame) => true
		@fact df[:a] => [0.0]
		@fact df[:b] => [0.0]

	end

	

end



end # module 






