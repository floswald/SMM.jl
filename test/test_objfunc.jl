module Test_objfunc

	using MOpt, Base.Test,DataStructures,DataFrames

	@testset "testing objfunctions" begin

		@testset "testing Testobj2" begin

			# pb    = OrderedDict("p1" => [0.2,-3,3] , "p2" => [-0.2,-2,2] )
			moms = DataFrame(name=["mu1","mu2"],value=[-1.0,1.0],weight=ones(2))
			mprob = MProb() 
			# addSampledParam!(mprob,pb) 
			addMoment!(mprob,moms) 
			# addEvalFunc!(mprob,objfunc_norm)

			ev = MOpt.Eval( mprob, OrderedDict(:p1 => 1.0 , :p2 => 0.0))
			@test ev.status == -1
			ev = MOpt.Testobj2(ev)	
			@test ev.status == 1
			@test param(ev,:p1) == 1.0
			@test param(ev,:p2) == 0.0


		end


		@testset "testing bivariate normal" begin

			ev = MOpt.Eval( Dict(:p1 => 1.0 , :p2 => 0.0), Dict(:mu1 =>0.0 , :mu2 => 0.0))
			ev = MOpt.objfunc_norm(ev)	
			@test abs(ev.simMoments[:mu1] - 1.0) < 0.1 
			@test abs(ev.simMoments[:mu2] - 1.0) > 0.1 

		end

	end


end # module 






