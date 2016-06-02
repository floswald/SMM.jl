module TestSlices

	using FactCheck,DataFrames,Lazy,MOpt

	facts("verifying read and write for Slices") do
		# initial value
		pb    = Dict("p1" => [0.0;-2;2] , "p2" => [0.0;-2;2] ) 
		moms = DataFrame(name=["mu1";"mu2"],value=[0.0;0.0],weight=rand(2))
		mprob = @> MOpt.MProb() MOpt.addSampledParam!(pb) MOpt.addMoment!(moms) MOpt.addEvalFunc!(MOpt.objfunc_norm);

		# compute the slices
		sl = MOpt.slices(mprob,30);

		MOpt.write(sl,"slices.h5")

		# get some slices for plotting
		r1 = MOpt.get(sl, :p1 , :mu1)
		r2 = MOpt.get(sl, :p2 , :value)

		# read from disk
		disk = MOpt.readSlice("slices.h5")
		# notice that i need to order the values here
		# this comes from floats being used as keys in the res. i store each param value as the (string) name of an HDF5 dataset, then convert it back to float.
		@fact all(sort(r1[:x]) .== sort(MOpt.get(disk,:p1,:mu1)[:x])) --> true
		@fact all(sort(r1[:y]) .== sort(MOpt.get(disk,:p1,:mu1)[:y])) --> true
		@fact all(sort(r2[:x]) .== sort(MOpt.get(disk,:p2,:value)[:x])) --> true
		@fact all(sort(r2[:y]) .== sort(MOpt.get(disk,:p2,:value)[:y])) --> true
	end
	rm("slices.h5")

	facts("handling of failing objective on remote") do

		context("can recover from RemoteException") do

			# include code here
			include(joinpath(Pkg.dir("MOpt"),"test","test-include.jl"))

			# add node and include code there.
			addprocs(1)
			remotecall_fetch(include,workers()[1],joinpath(Pkg.dir("MOpt"),"test","test-include.jl"))

			sl = MOpt.slices(mprob_fail,30);

			rmprocs(workers())

			# worked
			@fact true --> true

		end
	end
	# using PyPlot

	# subplot(231)
	# r = MOpt.get(sl, :p1 , :mu1); PyPlot.plot(r[:x],r[:y],".")
	# title("p1 vs mu1")
	# subplot(232)
	# r = MOpt.get(sl, :p1 , :mu2); PyPlot.plot(r[:x],r[:y],".")
	# title("p1 vs mu2")
	# subplot(233)
	# r = MOpt.get(sl, :p1 , :value); PyPlot.plot(r[:x],r[:y],".")
	# title("p1 vs value")
	# subplot(234)
	# r = MOpt.get(sl, :p2 , :mu1); PyPlot.plot(r[:x],r[:y],".")
	# title("p2 vs mu1")
	# subplot(235)
	# r = MOpt.get(sl, :p2 , :mu2); PyPlot.plot(r[:x],r[:y],".")
	# title("p2 vs mu2")
	# subplot(236)
	# r = MOpt.get(sl, :p2 , :value); PyPlot.plot(r[:x],r[:y],".")
	# title("p2 vs value")

end # module 






