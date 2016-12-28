module TestBGPChain

using Base.Test, DataFrames, MOpt, Lazy



# TESTING Chains
# ==============
pb   = Dict( "a" => [0.3; 0;1] , "b" => [4.9;0;1] )
moms = DataFrame(name=["alpha";"beta";"gamma"],value=[0.8;0.7;0.5],weight=rand(3))

@testset "Testing BGPChain constructor" begin
	
	mprob = MProb() 
	addSampledParam!(mprob,pb) 
	addMoment!(mprob,moms) 
	addEvalFunc!(mprob,MOpt.Testobj)
	id = 180
	n = 23
	sig = rand(length(pb))
	sig2 = rand(length(pb))
	upd = 5
	chain = MOpt.BGPChain(id,n,mprob,sig,upd)

	@test chain.id == id
	@test chain.accept_rate == 0.0
	@test chain.iter == 0
	@test chain.accepted == falses(n)
	@test chain.exchanged == zeros(Int,n)
	@test diag(chain.sigma) == sig
	@test chain.update_sigma == upd

	@testset "methods" begin
		chain.iter = 1
		ev = Eval(mprob)
		v = rand()
		ev.value = v
		ev.accepted = true
		MOpt.set_eval!(chain,ev)
		@test isa(chain.evals[1],Eval)
		@test chain.evals[1].value == v 
		@test chain.accepted[1]
		MOpt.set_acceptRate!(chain)
		@test chain.accept_rate == 1.0

		MOpt.set_sigma!(chain,sig2)
		@test diag(chain.sigma) == sig2

		MOpt.getLastEval(chain) == ev


	end


end



end # module 