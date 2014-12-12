module Test_objfunc

	using MOpt, FactCheck,DataFrames,Lazy

	facts("testing objfunctions") do

		ev = MOpt.Eval( [:m1 => 1.0 , :m2 => 0.0], [:m1 =>0.0 , :m2 => 0.0])
		ev = MOpt.objfunc_norm(ev)
		@fact abs(ev.moments[:m1] - 1.0) < 0.1 => true
		@fact abs(ev.moments[:m2] - 1.0) > 0.1 => true

		ev = MOpt.Eval( [:m1 => 0.0 , :m2 => 1.0], [:m1 =>0.0 , :m2 => 0.0])
		ev = MOpt.objfunc_norm(ev)
		@fact abs(ev.moments[:m1] - 1.0) > 0.1 => true
		@fact abs(ev.moments[:m2] - 1.0) < 0.1 => true

	end


end # module 






