
# TESTING Chains
# ==============


@testset "Testing BGPChains constructor" begin

	@testset "constructor" begin

	    (chain, id, n, mprob, sig,sig2,  upd, upd_by, ite, bundle) = test_chain()
		@test chain.id == id
		@test chain.accept_rate == 0.0
		@test chain.iter == 0
		@test chain.smpl_iters == ite
		@test chain.accepted == falses(n)
		@test chain.exchanged == zeros(Int,n)
		for (six,ix) in enumerate(chain.batches)
			@test chain.sigma == sig
		end
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
		SMM.set_eval!(chain,ev)
		@test isa(chain.evals[1],Eval)
		@test chain.evals[1].value == v 
		@test chain.accepted[1]
		SMM.set_acceptRate!(chain)
		@test chain.accept_rate == 1.0

		# SMM.set_sigma!(chain,sig2)
		# @test diag(chain.sigma) == sig2

		SMM.getLastAccepted(chain) == ev
	end

	@testset "proposal" begin
	    (chain, id, n, mprob, sig, sig2, upd, upd_by, ite) = test_chain()
		@test chain.iter == 0
		chain.iter = 1
		pp = SMM.proposal(chain)
		@test pp == mprob.initial_value

		#  set a value:
		ev = Eval(mprob)
		v = rand()
		ev.value = v
		ev.accepted = true
		SMM.set_eval!(chain,ev)
		# next period:
		chain.iter += 1
		pp = SMM.proposal(chain)
		@test pp != mprob.initial_value

	end

	# https://github.com/floswald/SMM.jl/issues/31
	@testset "high dimensional proposal" begin
	    (chain, id, n, mprob, sig, sig2, upd, upd_by, ite) = test_chain3()
		@test chain.iter == 0
		chain.iter = 1
		pp = SMM.proposal(chain)
		@test pp == mprob.initial_value

		#  set a value:
		ev = Eval(mprob)
		v = rand()
		ev.value = v
		ev.accepted = true
		# test whether we can get proposals out of this hidim MVNormal
		for i in 1:n-1
			ev.value = rand()
			ev.accepted = true
			SMM.set_eval!(chain,ev)
			chain.iter += 1
			p = SMM.proposal(chain)
			@test p != pp
			pp = deepcopy(p)
		end

	end

	@testset "testing accept reject" begin
	    (c, id, n, mprob, sig, sig2, upd, upd_by, ite) = test_chain()

		@testset "testing initial period" begin
			# nothing gets rejected in period 1
		    c.iter += 1

			# get a getNewCandidates
		    pp = SMM.proposal(c)
		    # evaluate objective 
		    ev = SMM.evaluateObjective(c.m,pp)

		    SMM.doAcceptReject!(c,ev)

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

			ev_0 = SMM.getLastAccepted(c)
			ev_good = Eval()
			ev_bad  = Eval()

			ev_good.value = 1.0
			ev_bad.value  = 2.0
			ev_good.status = 1
			ev_bad.status  = 1


			SMM.doAcceptReject!(c,ev_bad)
			@test ev_bad.prob < 1
			if ev_bad.prob > c.probs_acc[c.iter]
				@test ev_bad.accepted
				@test c.accepted[c.iter]
			end

			# force accept 
			ev_bad.accepted = true
			SMM.set_eval!(c,ev_bad)

			SMM.doAcceptReject!(c,ev_good)
			@test ev_good.prob == 1.0
			@test ev_good.accepted
			@test c.accepted[c.iter]

		end

	end
end



