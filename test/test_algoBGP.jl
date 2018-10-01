



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
	    o = MomentOpt.serialNormal(2,20)[1];
	    h = MomentOpt.history(o.chains[1])
	    @test isa(o,MAlgoBGP)
	    @test isa(h,DataFrame)
	    @test nrow(h) == 20
	    @test ncol(h) == 9
	    @test o.i == 20
	end

	@testset "parallelNormal() runs" begin
		# addprocs()

		addprocs(1)

		@everywhere using MomentOpt
	    o = MomentOpt.parallelNormal(20);
	    h = MomentOpt.history(o.chains[1])
	    @test isa(o,MAlgoBGP)
	    @test isa(h,DataFrame)
	    @test nrow(h) == 20
	    @test ncol(h) == 9
	    @test o.i == 20

	    rmprocs(workers())
	end

	@testset "Can we recover the mean of a Normal?" begin

		# set tolerance level
		# to achieve a smaller tolerance level, increase niter below
		# but takes more time
		tolTestNormal = 0.7

		addprocs(2)

		@everywhere using MomentOpt

		srand(1234)
		pb = OrderedDict("p1" => [0.2,-3,3] , "p2" => [-0.2,-2,2])
		trueValues = OrderedDict("mu1" => [-1.0] , "mu2" => [1.0])
		moms = DataFrame(name=["mu1","mu2"],value=[-1.0, 1.0], weight=ones(2))

		mprob = MProb()
		addSampledParam!(mprob,pb)
		addMoment!(mprob,moms)
		addEvalFunc!(mprob, MomentOpt.objfunc_norm)


		# estimation options:
		#--------------------
		niter = 200
		nchains = 2

		opts = Dict("N"=>nchains,
		        "maxiter"=>niter,
		        "maxtemp"=> 5,
		        "coverage"=>0.025,
		        "sigma_update_steps"=>10,
		        "sigma_adjust_by"=>0.01,
		        "smpl_iters"=>1000,
		        "parallel"=>true,
		        "maxdists"=>[0.0 for i in 1:nchains],
		        "acc_tuners"=>[5;1.0],
		        "animate"=>false)
		MA = MAlgoBGP(mprob,opts)
		MomentOpt.runMOpt!(MA)

		dat_chain1 = MomentOpt.history(MA.chains[1])
		dat_chain1[round(Int, niter/10):niter, :]
		dat_chain1 = dat_chain1[dat_chain1[:accepted ].== true, : ]

		parameters = [Symbol(String("mu$(i)")) for i=1:2]
		estimatedParameters = [Symbol(String("p$(i)")) for i=1:2]

		for (estimatedParameter, param) in zip(estimatedParameters, parameters)

		  println("Quasi posterior mean for $(String(estimatedParameter)) = $(mean(dat_chain1[estimatedParameter]))")
		  println("Quasi posterior median for $(String(estimatedParameter)) = $(median(dat_chain1[estimatedParameter]))")
		  println("True value = $(trueValues[String(param)][])")

			# (the problem is simple here)
			#---------------------------------------------------------
			@test trueValues[String(param)][] ≈ median(dat_chain1[estimatedParameter]) atol = tolTestNormal

		end

	    rmprocs(workers())
	end

	@testset "with batches?" begin

		# set tolerance level
		# to achieve a smaller tolerance level, increase niter below
		# but takes more time
		tolTestNormal = 0.7

		addprocs(2)

		@everywhere using MomentOpt

		srand(1234)
		pb = OrderedDict("p1" => [0.2,-3,3] , "p2" => [-0.2,-2,2])
		trueValues = OrderedDict("mu1" => [-1.0] , "mu2" => [1.0])
		moms = DataFrame(name=["mu1","mu2"],value=[-1.0, 1.0], weight=ones(2))

		mprob = MProb()
		addSampledParam!(mprob,pb)
		addMoment!(mprob,moms)
		addEvalFunc!(mprob, MomentOpt.objfunc_norm)


		# estimation options:
		#--------------------
		niter = 50
		nchains = 2

		opts = Dict("N"=>nchains,
		        "maxiter"=>niter,
		        "maxtemp"=> 5,
		        "coverage"=>0.025,
		        "sigma_update_steps"=>10,
		        "sigma_adjust_by"=>0.01,
		        "smpl_iters"=>1000,
		        "parallel"=>true,
		        "maxdists"=>[0.05 for i in 1:nchains],
		        "mixprob"=>0.3,
		        "acc_tuner"=>12.0,
		        "animate"=>false,
		        "batch_size"=> 1)
		MA = MAlgoBGP(mprob,opts)
		MomentOpt.runMOpt!(MA)

		dat_chain1 = MomentOpt.history(MA.chains[1])
		dat_chain1[round(Int, niter/10):niter, :]
		dat_chain1 = dat_chain1[dat_chain1[:accepted ].== true, : ]

		parameters = [Symbol(String("mu$(i)")) for i=1:2]
		estimatedParameters = [Symbol(String("p$(i)")) for i=1:2]

		for (estimatedParameter, param) in zip(estimatedParameters, parameters)

		  println("Quasi posterior mean for $(String(estimatedParameter)) = $(mean(dat_chain1[estimatedParameter]))")
		  println("Quasi posterior median for $(String(estimatedParameter)) = $(median(dat_chain1[estimatedParameter]))")
		  println("True value = $(trueValues[String(param)][])")

			# (the problem is simple here)
			#---------------------------------------------------------
			@test trueValues[String(param)][] ≈ median(dat_chain1[estimatedParameter]) atol = tolTestNormal

		end

	    rmprocs(workers())
	end

	# Running the estimation for 10 iterations and then adding 10 more iterations
	# should give the same results as running the estimation for 20 iterations
	# in the first place (keeping random shocks constant)
	@testset "Testing stopping and restarting" begin

	
		pb = OrderedDict("p1" => [0.2,-3,3] , "p2" => [-0.2,-2,2])
		trueValues = OrderedDict("mu1" => [-1.0] , "mu2" => [1.0])
		moms = DataFrame(name=["mu1","mu2"],value=[-1.0, 1.0], weight=ones(2))

		
		mprob = MProb()
		addSampledParam!(mprob,pb)
		addMoment!(mprob,moms)

		addEvalFunc!(mprob, objfunc_norm)

		niter = 20
		nchains = 2

		opts = Dict("N"=>nchains,
		        "maxiter"=>niter,
		        "maxtemp"=> 5,
		        "coverage"=>0.025,
		        "sigma_update_steps"=>10,
		        "sigma_adjust_by"=>0.01,
		        "smpl_iters"=>1000,
		        "parallel"=>true,
		        "maxdists"=>[0.05 for i in 1:nchains],
		        "mixprob"=>0.3,
		        "acc_tuner"=>12.0,
		        "animate"=>false)


		# set-up BGP algorithm:
		srand(1234)
		MA = MAlgoBGP(mprob,opts)

		# run the estimation
		MomentOpt.runMOpt!(MA)

		# restart the estimation:
		newiters = 20
		MomentOpt.restartMOpt!(MA, newiters)

		mprob2 = MProb()
		addSampledParam!(mprob2,pb)
		addMoment!(mprob2,moms)
		addEvalFunc!(mprob2, objfunc_norm)
		niter = niter + 	newiters
		# number of chains
		nchains = 2

		opts2 = Dict("N"=>nchains,
		        "maxiter"=>niter,
		        "maxtemp"=> 5,
		        "coverage"=>0.025,
		        "sigma_update_steps"=>10,
		        "sigma_adjust_by"=>0.01,
		        "smpl_iters"=>1000,
		        "parallel"=>true,
		        "maxdists"=>[0.05 for i in 1:nchains],
		        "mixprob"=>0.3,
		        "acc_tuner"=>12.0,
		        "animate"=>false)


		# set-up BGP algorithm:
		srand(1234)
		MA2 = MAlgoBGP(mprob2,opts2)

		# run the estimation:
		MomentOpt.runMOpt!(MA2)

		MomentOpt.summary(MA2)
		MomentOpt.summary(MA)

		Qmean = zeros(2,2)
		Qmedian = zeros(2,2)

		for (indexAlgo, algo) in enumerate([MA, MA2])

			# Realization of the first chain:
			#-------------------------------
			dat_chain1 = MomentOpt.history(algo.chains[1])

			# discard the first 10th of the observations ("burn-in" phase):
			#--------------------------------------------------------------
			dat_chain1[round(Int, niter/10):niter, :]

			# keep only accepted draws:
			#--------------------------
			dat_chain1 = dat_chain1[dat_chain1[:accepted ].== true, : ]

			# create a list with the parameters to be estimated
			parameters = [Symbol(String("mu$(i)")) for i=1:2]
			# list with the corresponding priors:
			#------------------------------------
			estimatedParameters = [Symbol(String("p$(i)")) for i=1:2]

			# Quasi Posterior mean and quasi posterior median for each parameter:
			#-------------------------------------------------------------------
			for (estimatedParameter, param) in zip(estimatedParameters, parameters)

			  Qmean[indexAlgo,:] = mean(dat_chain1[estimatedParameter])
			  Qmedian[indexAlgo,:] = median(dat_chain1[estimatedParameter])

			end

		end

		# Compare values
		#---------------
		@test Qmean[1,1] .== Qmean[1,2]
		@test Qmean[2,1] .== Qmean[2,2]
		@test Qmedian[1,1] .== Qmedian[1,2]
		@test Qmedian[2,1] .== Qmedian[2,2]

	end


end
