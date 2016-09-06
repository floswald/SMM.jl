module Test_objfunc

	using MOpt, FactCheck,DataFrames,Lazy

	facts("testing objfunctions") do

		context("testing Testobj2") do

			ev = MOpt.Eval( Dict(:p1 => 1.0 , :p2 => 0.0), Dict(:mu1 =>0.0 , :mu2 => 0.0))
			@fact ev.status --> -1
			ev = MOpt.Testobj2(ev)	
			@fact ev.status --> 1
			@fact param(ev,:p1) --> 1.0
			@fact param(ev,:p2) --> 0.0


		end


		context("testing bivariate normal") do

			ev = MOpt.Eval( Dict(:p1 => 1.0 , :p2 => 0.0), Dict(:mu1 =>0.0 , :mu2 => 0.0))
			ev = MOpt.objfunc_norm(ev)	
			@fact abs(ev.simMoments[:mu1] - 1.0) < 0.1 --> true
			@fact abs(ev.simMoments[:mu2] - 1.0) > 0.1 --> true

		end

	end


end # module 






