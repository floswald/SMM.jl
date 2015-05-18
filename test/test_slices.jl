module TestSlices

	using FactCheck,DataFrames,Lazy,MOpt

	# initial value
	pb    = ["m1" => [0.0,-2,2] , "m2" => [0.0,-2,2] ] 
	moms = DataFrame(name=["mu1","mu2"],value=[0.0,0.0],weight=rand(2))
	mprob = @> MOpt.MProb() MOpt.addSampledParam!(pb) MOpt.addMoment!(moms) MOpt.addEvalFunc!(MOpt.objfunc_norm);

	# compute the slices
	sl = MOpt.slices(mprob,30);

	# using PyPlot

	# subplot(231)
	# r = MOpt.get(sl, :m1 , :m1); PyPlot.plot(r[:x],r[:y],".")
	# subplot(232)
	# r = MOpt.get(sl, :m1 , :m2); PyPlot.plot(r[:x],r[:y],".")
	# subplot(233)
	# r = MOpt.get(sl, :m1 , :value); PyPlot.plot(r[:x],r[:y],".")
	# subplot(234)
	# r = MOpt.get(sl, :m2 , :m1); PyPlot.plot(r[:x],r[:y],".")
	# subplot(235)
	# r = MOpt.get(sl, :m2 , :m2); PyPlot.plot(r[:x],r[:y],".")
	# subplot(236)
	# r = MOpt.get(sl, :m2 , :value); PyPlot.plot(r[:x],r[:y],".")

end # module 






