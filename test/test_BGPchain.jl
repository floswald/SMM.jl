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
		upd_by = rand()
		ite = 10
		chain = MOpt.BGPChain(id,n,mprob,sig,upd,upd_by,ite)

		@test chain.id == id
		@test chain.accept_rate == 0.0
		@test chain.iter == 0
		@test chain.smpl_iters == ite
		@test chain.accepted == falses(n)
		@test chain.exchanged == zeros(Int,n)
		@test diag(chain.sigma) == sig
		@test chain.sigma_update_steps == upd
		@test chain.sigma_adjust_by == upd_by

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
		upd_by = rand()
		ite = 10
		chain = MOpt.BGPChain(id,n,mprob,sig,upd,upd_by,ite)
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
		upd_by = rand()
		ite = 5000
		chain = MOpt.BGPChain(id,n,mprob,sig,upd,upd_by,ite)
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

	@testset "testing accept reject" begin


		@testset "testing initial period" begin
			# nothing gets rejected in period 1
			c = MOpt.BGPChain()

		    # increment interation
		    c.iter += 1

			# get a getNewCandidates
		    pp = MOpt.getNewCandidates(c)
		    # evaluate objective 
		    ev = MOpt.evaluateObjective(c.m,pp)

		    MOpt.doAcceptReject!(c,ev)

			# is accepted: 
			@test c.accepted[c.iter]
			@test ev.prob == 1
			@test ev.accepted

		end

		# @testset "testing whether params get accept/rejected" begin

		# 	# first iteration
		# 	MA = MAlgoBGP(mprob,opts)
		# 	MA.i = 1
		# 	MOpt.incrementChainIter!(MA.MChains)
		# 	v = map( x -> MOpt.evaluateObjective(MA.m,x), MA.current_param)

		# 	MOpt.doAcceptRecject!(MA,v)

		# 	# second iteration
		# 	MA.i = 2
		# 	MOpt.incrementChainIter!(MA.MChains)

		# 	# get randome parameters
		# 	for i in 1:length(v) 
		# 		for (k,va) in MA.current_param[i]
		# 			MA.current_param[i][k] = va + rand()
		# 		end
		# 	end

		# 	# get 2 return values
		# 	v1 = map( x -> MOpt.evaluateObjective(MA.m,x), MA.current_param)
		# 	v2 = map( x -> MOpt.evaluateObjective(MA.m,x), MA.current_param)

		# 	# assign a very good bad value on each chain: 
		# 	for i in 1:length(v) 
		# 		v1[i].value = v1[i].value - 100.0
		# 		setMoment(v1[i],  Dict( :alpha => rand(), :beta => rand(), :gamma => rand() ))
		# 		v2[i].value = v2[i].value + 100.0
		# 		setMoment(v2[i],  Dict( :alpha => rand(), :beta => rand(), :gamma => rand() ) )
		# 	end
		# 	MAs = MAlgo[]
		# 	push!(MAs,deepcopy(MA),deepcopy(MA))
		# 	MOpt.doAcceptRecject!(MAs[1],v1)
		# 	MOpt.doAcceptRecject!(MAs[1],v2)

		# 	for ma in MAs
		# 		for ch in 1:ma["N"]
		# 			ev= getEval(ma.MChains[ch],ma.i)
		# 			Last_ev= getEval(ma.MChains[ch],ma.i-1)

		# 			if infos(ma.MChains[ch],ma.i)[:accept][1]
		# 				# if accepted, the parameter vector in ev is equal to current_param on that chain
		# 				ev_p = paramd(ev)
		# 				ev_m = ev.simMoments
		# 				for (k,v) in ev_p
		# 					@test v == ma.current_param[ch][k]
		# 				end
		# 			else
		# 				# if not, parameters(ma.MChains[ch],ma.i) == parameters(ma.MChains[ch],ma.i-1)
		# 				ev_p = paramd(Last_ev)
		# 				ev_m = Last_ev.simMoments
		# 				for (k,v) in ev_p
		# 					@test v == not(ma.current_param[ch][k])
		# 				end
		# 			end
		# 		end
		# 	end
		# end
	end


























end



end # module 