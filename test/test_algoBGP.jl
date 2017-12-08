



@testset "AlgoBGP" begin


	pb   = Dict( "a" => [0.3; -1;1] , "b" => [-0.9;-2;2] )
	moms = DataFrame(name=["mu1","mu2"],value=[0.0;0.0],weight=rand(2))
	mprob = MProb() 
	addSampledParam!(mprob,pb) 
	addMoment!(mprob,moms) 
	addEvalFunc!(mprob,MomentOpt.objfunc_norm)
	@testset "Constructor" begin

		MA = MAlgoBGP(mprob)

		@test isa(MA,MAlgo) == true
		@test isa(MA,MAlgoBGP) == true
		@test isa(MA.m,MProb) == true

		@test MA.i == 0
		@test length(MA.chains) == 3

		for ix = 1:length(MA.chains)
			@test isa( MA.chains[ix] , MomentOpt.BGPChain)
		end
	end

	@testset "serialNormal() runs" begin
	    o = MomentOpt.serialNormal(20);
	    h = MomentOpt.history(o.chains[1])
	    @test isa(o,MAlgoBGP)
	    @test isa(h,DataFrame)
	    @test nrow(h) == 20
	    @test ncol(h) == 9
	    @test o.i == 20
	end

	
end


