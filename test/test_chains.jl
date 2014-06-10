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
		@fact length(chain.infos["evals"])  => L
		@fact length(chain.infos["accept"]) => L
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

facts("testing MChain constructor") do
	

	mprob = Mopt.MProb(p,pb,Mopt.Testobj,moms)
	L = 9
	chain = Mopt.Chain(mprob,L)

	n = 14
	chains = Mopt.MChain(n,mprob,L)

	@fact isa(chains,Mopt.MChain) => true
	@fact isa(chains.chains,Array) => true
	@fact length(chains.chains) => n
	for ch in 1:n
		@fact length(chains.chains[ch].infos["evals"]) => L
	end
end

facts("testing Chain/MChain methods") do

	context("testing getindex(chain,i)") do

		mprob = Mopt.MProb(p,pb,Mopt.Testobj,moms)
		v = Mopt.Testobj(p,moms,["alpha","beta","gamma"])
		L = 9	# length of chain
		chain = Mopt.Chain(mprob,L)

		i = rand(1:L)

		x = Mopt.alls(chain,i,true)

		@fact isa(x,DataFrame) => true
		@fact nrow(x) => 1
		@fact ncol(x) => 2 + length(Mopt.ps_names(mprob)) + length(Mopt.ms_names(mprob))

		nx = names(x)
		nm = [:evals,:accept,map(x-> symbol(x), collect(Mopt.ps_names(mprob))),map(x-> symbol(x), collect(Mopt.ms_names(mprob)))]
		@fact sort(nx) == sort(nm) => true

	end

	context("testing getindex(chain,i::Range)") do

		mprob = Mopt.MProb(p,pb,Mopt.Testobj,moms)
		v     = Mopt.Testobj(p,moms,["alpha","beta","gamma"])
		L     = 9	# length of chain
		chain = Mopt.Chain(mprob,L)

		i = 3:L
		x = Mopt.alls(chain,i,true)
		
		@fact isa(x,DataFrame) => true
		@fact nrow(x) => length(i)
		@fact ncol(x) => 2 + length(Mopt.ps_names(mprob)) + length(Mopt.ms_names(mprob))

		nx = names(x)
		nm = [:evals,:accept,map(x-> symbol(x), collect(Mopt.ps_names(mprob))),map(x-> symbol(x), collect(Mopt.ms_names(mprob)))]
		@fact sort(nx) == sort(nm) => true

	end
	
	context("test appendEval!(chain)") do

		mprob = Mopt.MProb(p,pb,Mopt.Testobj,moms)
		v = Mopt.Testobj(p,moms,["alpha","beta","gamma"])
		L = 9
		chain = Mopt.Chain(mprob,L)

		# verify values are zero:
		@fact all(chain.infos["evals"] .== 0.0) => true
		@fact chain.i => 0 

		# set i to 1 to test this:
		chain.i = 1

		# update chain with v
		Mopt.appendEval!(chain,v,true)

		# verify new values on chain
		@fact chain.i => 1 
		@fact chain.infos["evals"][1] => v["value"]
		@fact chain.infos["accept"][1] => true
		for nm in Mopt.ps_names(mprob)
			@fact chain.parameters[nm][1] => v["params"][nm]
		end
		for nm in Mopt.ms_names(mprob)
			@fact chain.moments[nm][1] => v["moments"][nm]
		end
	end

	context("test appendEval!(chains)") do

		mprob = Mopt.MProb(p,pb,Mopt.Testobj,moms)
		L = 9
		n = 17
		v = [ Mopt.Testobj(p,moms,["alpha","beta","gamma"]) for i=1:n ]
		MC = Mopt.MChain(n,mprob,L)

		# verify values are zero:
		@fact all(MC.chains[1].infos["evals"] .== 0.0) => true
		@fact all(MC.chains[1].infos["evals"] .== 0.0) => true
		@fact MC.chains[1].i => 0 

		@fact isempty(Mopt.evals(MC[1])) => true

		# set i to 1 to test this, otherwise BoundsError (accessing x[0])
		for ix = 1:n
			MC[ix].i = 1
		end

		# update chain with v
		which = rand(1:n)
		Mopt.appendEval!(MC,which,v[which],true)

		# verify new values on each chain
		@fact MC[which].i => 1 
		@fact Mopt.evals(MC[which])[1] => v[which]["value"]
		@fact MC.chains[which].infos["accept"][1] => true
		for nm in Mopt.ps_names(mprob)
			@fact Mopt.parameters(MC[which])[nm][1] => v[which]["params"][nm]
		end
		for nm in Mopt.ms_names(mprob)
			@fact Mopt.moments(MC[which])[nm][1] => v[which]["moments"][nm]
		end
	end

	context("testing updateIter(MChain") do

		mprob = Mopt.MProb(p,pb,Mopt.Testobj,moms)
		L = 9
		n = 17
		MC = Mopt.MChain(n,mprob,L)

		for ix = 1:n
			@fact MC.chains[ix].i => 0
		end

		Mopt.updateIter!(MC)
		for ix = 1:n
			@fact MC.chains[ix].i => 1
		end

	end


end





end # module 






