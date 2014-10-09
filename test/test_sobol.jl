module TestSobol

	using FactCheck,DataFrames,Lazy,MOpt

	# initial value
	pb    = ["m1" => [0.0,-2,2] , "m2" => [0.0,-2,2] ] 
	moms = DataFrame(name=["m1","m2"],value=[0.2,-0.6],weight=rand(2))
	mprob = @> MOpt.MProb() MOpt.addSampledParam!(pb) MOpt.addMoment(moms) MOpt.addEvalFunc(MOpt.objfunc_norm);

	# compute the slices
	sl = MOpt.sobolsearch(mprob,30);

	# find the min value

	best=sl[1]
	for e in sl
		if e.value < best.value
			best=e
		end
	end

	print(best.params)

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






