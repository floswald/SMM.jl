
pb    = Dict("p1" => [0.0;-2;2] , "p2" => [0.0;-2;2] ) 
moms = DataFrame(name=["mu1";"mu2"],value=[0.0;0.0],weight=rand(2))
mprob = MomentOpt.MProb() 
MomentOpt.addSampledParam!(mprob,pb) 
MomentOpt.addMoment!(mprob,moms) 
MomentOpt.addEvalFunc!(mprob,MomentOpt.objfunc_norm)
mprob_fail = deepcopy(mprob)
mprob_fail.objfunc = MomentOpt.Testobj_fails

@testset "verifying read and write for Slices" begin
	# initial value

	# compute the slices
	sl = MomentOpt.doSlices(mprob,30);

	MomentOpt.write(sl,"slices.h5")

	# get some slices for plotting
	r1 = MomentOpt.get(sl, :p1 , :mu1)
	r2 = MomentOpt.get(sl, :p2 , :value)

	# read from disk
	disk = MomentOpt.readSlice("slices.h5")
	# notice that i need to order the values here
	# this comes from floats being used as keys in the res. i store each param value as the (string) name of an HDF5 dataset, then convert it back to float.
	@test all(sort(r1[:x]) .== sort(MomentOpt.get(disk,:p1,:mu1)[:x])) 
	@test all(sort(r1[:y]) .== sort(MomentOpt.get(disk,:p1,:mu1)[:y])) 
	@test all(sort(r2[:x]) .== sort(MomentOpt.get(disk,:p2,:value)[:x]))
	@test all(sort(r2[:y]) .== sort(MomentOpt.get(disk,:p2,:value)[:y]))
	if !is_windows() rm("slices.h5") end
end

@testset "handling of failing objective on remote" begin

	@testset "can recover from RemoteException" begin

		# include code here
		# include(joinpath(Pkg.dir("MOpt"),"test","test-include.jl"))

		# add node and include code there.
		addprocs(1)
		remotecall_fetch(include,workers()[1],joinpath(dirname(@__FILE__),"..","src","MomentOpt.jl"))

		sl = MomentOpt.doSlices(mprob_fail,30);

		rmprocs(workers())

		# worked
		@test true

	end
end
# using PyPlot

# subplot(231)
# r = MomentOpt.get(sl, :p1 , :mu1); PyPlot.plot(r[:x],r[:y],".")
# title("p1 vs mu1")
# subplot(232)
# r = MomentOpt.get(sl, :p1 , :mu2); PyPlot.plot(r[:x],r[:y],".")
# title("p1 vs mu2")
# subplot(233)
# r = MomentOpt.get(sl, :p1 , :value); PyPlot.plot(r[:x],r[:y],".")
# title("p1 vs value")
# subplot(234)
# r = MomentOpt.get(sl, :p2 , :mu1); PyPlot.plot(r[:x],r[:y],".")
# title("p2 vs mu1")
# subplot(235)
# r = MomentOpt.get(sl, :p2 , :mu2); PyPlot.plot(r[:x],r[:y],".")
# title("p2 vs mu2")
# subplot(236)
# r = MomentOpt.get(sl, :p2 , :value); PyPlot.plot(r[:x],r[:y],".")
# title("p2 vs value")







