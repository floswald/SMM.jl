module TestBGPChain

using Base.Test, DataFrames, MOpt, Lazy



# TESTING Chains
# ==============
pb   = Dict( "a" => [0.3; 0;1] , "b" => [0.9;0;1] )
moms = DataFrame(name=["alpha";"beta";"gamma"],value=[0.8;0.7;0.5],weight=rand(3))

@testset "Testing BGPChains constructor" begin
	
	@testset "constructor" begin
		mprob = MProb() 
		addSampledParam!(mprob,pb) 
		addMoment!(mprob,moms) 
		addEvalFunc!(mprob,MOpt.Testobj)
		id = 180
		n = 23
		sig = rand(length(pb))
		sig2 = rand(length(pb))
		upd = 5
		ite = 10
		chain = MOpt.BGPChain(id,n,mprob,sig,upd,ite)

		@test chain.id == id
		@test chain.accept_rate == 0.0
		@test chain.iter == 0
		@test chain.smpl_iters == ite
		@test chain.accepted == falses(n)
		@test chain.exchanged == zeros(Int,n)
		@test diag(chain.sigma) == sig
		@test chain.update_sigma == upd

	end

	@testset "basic methods" begin
		mprob = MProb() 
		addSampledParam!(mprob,pb) 
		addMoment!(mprob,moms) 
		addEvalFunc!(mprob,MOpt.Testobj)
		id = 180
		n = 23
		sig = rand(length(pb))
		sig2 = rand(length(pb))
		upd = 5
		ite = 10
		chain = MOpt.BGPChain(id,n,mprob,sig,upd,ite)
		@test chain.iter == 0
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

	@testset "sample similar variances" begin
		n = 10
		sig = rand(n)
		d = MOpt.MvNormal(zeros(10),MOpt.PDiagMat(sig))
		lb = [-0.5 for i in 1:n]
		ub = [ 0.5 for i in 1:n]
		x = MOpt.sample(d,lb,ub,1000)
		@test length(x)==n
 	end
	@testset "sample non-similar variances" begin
		n = 10
		sig = collect(linspace(1.0,10.0,n))
		d = MOpt.MvNormal(zeros(10),MOpt.PDiagMat(sig))
		lb = -2 * sig
		ub =  2 * sig
		x = MOpt.sample(d,lb,ub,1000)
		@test length(x)==n
 	end

	@testset "getNewCandidates" begin
		mprob = MProb() 
		addSampledParam!(mprob,pb) 
		addMoment!(mprob,moms) 
		addEvalFunc!(mprob,MOpt.Testobj)
		id = 180
		n = 23
		sig = rand(length(pb))
		upd = 5
		ite = 5000
		chain = MOpt.BGPChain(id,n,mprob,sig,upd,ite)
		@test chain.iter == 0
		chain.iter = 1
		pp = MOpt.getNewCandidates(chain)
		@test pp == mprob.initial_value

		#  set a value:
		ev = Eval(mprob)
		v = rand()
		ev.value = v
		ev.accepted = true
		MOpt.set_eval!(chain,ev)
		# next period:
		chain.iter += 1
		pp = MOpt.getNewCandidates(chain)
		@test pp != mprob.initial_value

	end



























end



end # module 