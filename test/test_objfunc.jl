@testset "testing objfunctions" begin

	@testset "testing Testobj2" begin

		# pb    = OrderedDict("p1" => [0.2,-3,3] , "p2" => [-0.2,-2,2] )
		moms = DataFrame(name=["mu1","mu2"],value=[-1.0,1.0],weight=ones(2))
		mprob = MProb() 
		# addSampledParam!(mprob,pb) 
		addMoment!(mprob,moms) 
		# addEvalFunc!(mprob,objfunc_norm)

		ev = MomentOpt.Eval( mprob, OrderedDict(:p1 => 1.0 , :p2 => 0.0))
		@test ev.status == -1
		ev = MomentOpt.Testobj2(ev)	
		@test ev.status == 1
		@test param(ev,:p1) == 1.0
		@test param(ev,:p2) == 0.0


	end

	@testset "testing bivariate normal" begin

		ev = MomentOpt.Eval( Dict(:p1 => 0.0 , :p2 => 0.0), Dict(:mu1 =>0.0 , :mu2 => 0.0))
		ev = MomentOpt.objfunc_norm(ev)	
		@test isapprox(ev.simMoments[:mu1], ev.dataMoments[:mu1], atol= 0.1 )
		@test isapprox(ev.simMoments[:mu2], ev.dataMoments[:mu2], atol= 0.1 )

	end
end








