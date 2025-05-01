
@testset "slices" begin
    

	pb    = Dict("p1" => [0.0;-2;2] , "p2" => [0.0;-2;2] ) 
	moms = DataFrame(name=["mu1";"mu2"],value=[0.0;0.0],weight=rand(2))
	mprob = SMM.MProb() 
	SMM.addSampledParam!(mprob,pb) 
	SMM.addMoment!(mprob,moms) 
	SMM.addEvalFunc!(mprob,SMM.objfunc_norm)
	mprob_fail = deepcopy(mprob)
	mprob_fail.objfunc = SMM.Testobj_fails

	@testset "verifying read and write for Slices" begin
		# initial value

		# compute the slices
		sl = SMM.doSlices(mprob,3);

		SMM.save(sl,"slices.jld2")

		# get some slices for plotting
		r1 = SMM.get(sl, :p1 , :mu1)
		r2 = SMM.get(sl, :p2 , :value)

		# read from disk
		disk = SMM.load("slices.jld2")["s"]
		# notice that i need to order the values here
		# this comes from floats being used as keys in the res. i store each param value as the (string) name of an HDF5 dataset, then convert it back to float.
		@test all(sort(r1[:x]) .== sort(SMM.get(disk,:p1,:mu1)[:x])) 
		@test all(sort(r1[:y]) .== sort(SMM.get(disk,:p1,:mu1)[:y])) 
		@test all(sort(r2[:x]) .== sort(SMM.get(disk,:p2,:value)[:x]))
		@test all(sort(r2[:y]) .== sort(SMM.get(disk,:p2,:value)[:y]))
		if !Sys.iswindows() rm("slices.jld2") end
	end

	@testset "handling of failing objective on remote" begin

		@testset "can recover from RemoteException" begin

			# include code here
			# include(joinpath(Pkg.dir("MOpt"),"test","test-include.jl"))

			# add node and include code there.
			addprocs(1)
			remotecall_fetch(include,workers()[1],joinpath(dirname(@__FILE__),"..","src","SMM.jl"))

			sl = SMM.doSlices(mprob_fail,30);

			rmprocs(workers())

			# worked
			@test true

		end
	end

	@testset "naive coordinate descent works" begin
	    m,s = SMM.snorm_6_taxi(0.1)
	    @test norm(m[!,:value] .- collect(values(s[:best][:p]))) < 0.5
	end

	@testset "naive coordinate descent works in parallel" begin
		addprocs(2)
		@everywhere using SMM
	    m,s = SMM.snorm_6_taxi(0.1,par=true)
		rmprocs(workers())
	    @test norm(m[!,:value] .- collect(values(s[:best][:p]))) < 0.5
	end

	@testset "simple slice optimization works" begin
	    x = SMM.sliceOpt(0.001)
	    @test abs(x[:optimized][:best][:p][:p1] - x[:targets][1,2]) < 0.2 
	    @test abs(x[:optimized][:best][:p][:p2] - x[:targets][2,2]) < 0.2 

	end






end
