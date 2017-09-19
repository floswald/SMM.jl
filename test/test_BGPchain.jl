
# TESTING Chains
# ==============


@testset "Testing BGPChains constructor" begin

	@testset "constructor" begin

	    (chain, id, n, mprob, sig,sig2,  upd, upd_by, ite) = test_chain()
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
	    (chain, id, n, mprob, sig, sig2, upd, upd_by, ite) = test_chain()
		@test chain.iter == 0
		chain.iter = 1
		ev = Eval(mprob)
		v = rand()
		ev.value = v
		ev.accepted = true
		MomentOpt.set_eval!(chain,ev)
		@test isa(chain.evals[1],Eval)
		@test chain.evals[1].value == v 
		@test chain.accepted[1]
		MomentOpt.set_acceptRate!(chain)
		@test chain.accept_rate == 1.0

		MomentOpt.set_sigma!(chain,sig2)
		@test diag(chain.sigma) == sig2

		MomentOpt.getLastAccepted(chain) == ev
	end

	# @testset "sample similar variances" begin
	# 	n = 10
	# 	sig = rand(n)
	# 	d = MomentOpt.MvNormal(zeros(10),MomentOpt.PDiagMat(sig))
	# 	lb = [-0.5 for i in 1:n]
	# 	ub = [ 0.5 for i in 1:n]
	# 	x = MomentOpt.mysample(d,lb,ub,1000)
	# 	@test length(x)==n
 # 	end
	# @testset "sample non-similar variances" begin
	# 	n = 10
	# 	sig = collect(linspace(1.0,10.0,n))
	# 	d = MomentOpt.MvNormal(zeros(10),MomentOpt.PDiagMat(sig))
	# 	lb = -2 * sig
	# 	ub =  2 * sig
	# 	x = MomentOpt.mysample(d,lb,ub,1000)
	# 	@test length(x)==n
 # 	end

	@testset "proposal" begin
	    (chain, id, n, mprob, sig, sig2, upd, upd_by, ite) = test_chain()
		@test chain.iter == 0
		chain.iter = 1
		pp = MomentOpt.proposal(chain)
		@test pp == mprob.initial_value

		#  set a value:
		ev = Eval(mprob)
		v = rand()
		ev.value = v
		ev.accepted = true
		MomentOpt.set_eval!(chain,ev)
		# next period:
		chain.iter += 1
		pp = MomentOpt.proposal(chain)
		@test pp != mprob.initial_value

	end

	@testset "testing accept reject" begin
	    (c, id, n, mprob, sig, sig2, upd, upd_by, ite) = test_chain()

		@testset "testing initial period" begin
			# nothing gets rejected in period 1
		    c.iter += 1

			# get a getNewCandidates
		    pp = MomentOpt.proposal(c)
		    # evaluate objective 
		    ev = MomentOpt.evaluateObjective(c.m,pp)

		    MomentOpt.doAcceptReject!(c,ev)

			# is accepted: 
			@test c.accepted[c.iter]
			@test ev.prob == 1
			@test ev.accepted

		end

		@testset "test correct accept/reject" begin

		    c.iter += 1
		    @test c.iter == 2
			# get 2 Evals: one with good, one with bad value
			# want to accept good and reject bad.

			ev_0 = MomentOpt.getLastAccepted(c)
			ev_good = Eval()
			ev_bad  = Eval()

			ev_good.value = ev_0.value - 10.0
			ev_bad.value  = ev_0.value + 10.0
			ev_good.status = 1
			ev_bad.status  = 1


			MomentOpt.doAcceptReject!(c,ev_bad)
			@test ev_bad.prob < 1
			if ev_bad.prob > c.probs_acc[c.iter]
				@test ev_bad.accepted
				@test c.accepted[c.iter]
			end

			MomentOpt.doAcceptReject!(c,ev_good)
			@test ev_good.prob == 1.0
			@test ev_good.accepted
			@test c.accepted[c.iter]

		end

	end
end



