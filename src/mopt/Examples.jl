	


function serialNormal(logmode="debug")
	# data are generated from a bivariate normal
	# with mu = [a,b] = [0,0]
	# aim: 
	# 1) sample [a',b'] from a space [-1,1] x [-1,1] and
	# 2) find true [a,b] by computing distance(S([a',b']), S([a,b]))
	#    and accepting/rejecting [a',b'] according to BGP
	# 3) S([a,b]) returns a summary of features of the data

	# initial value
	pb    = Dict("p1" => [0.2,-2,2] , "p2" => [-0.2,-2,2] )
	moms = DataFrame(name=["mu2","mu1"],value=[-1.0,1.0],weight=ones(2))
	mprob = MProb() 
	addSampledParam!(mprob,pb) 
	addMoment!(mprob,moms) 
	addEvalFunc!(mprob,objfunc_norm)

	opts =Dict("N"=>3,
		"maxiter"=>100,
		"maxtemp"=> 2,
		"bound_prob"=>0.15,
		"disttol"=>0.00,
		"sigma_update_steps"=>1,
		"sigma_adjust_by"=>0.1,
		"smpl_iters"=>1000,
		"parallel"=>false)

	# setup the BGP algorithm
	MA = MAlgoBGP(mprob,opts)
	# MOpt.cur_param(MA)
	plot(MA,1)

	# run the estimation
	# runMOpt!(MA)
	# fig = figure("parameter histograms") 
	# plt[:hist](convert(Array,MOpt.parameters(MA.MChains[1])),15)
	# return MA
end


function BGP_example()

	p = Dict("theta" => [2.0,-2,10])
	mprob = MProb()
	addSampledParam!(mprob,p)
	addEvalFunc!(mprob,MOpt.objfunc_BGP)

	opts =Dict(
		"N"               => length(workers()),							# number of MCMC chains
		"maxiter"         => 500,						# max number of iterations
		"savefile"        => joinpath(pwd(),"MA.h5"),	# filename to save results
		"print_level"     => 1,							# increasing verbosity level of output
		"maxtemp"         => 1,							# tempering of hottest chain
		"min_shock_sd"    => 0.1,						# initial sd of shock on coldest chain
		"max_shock_sd"    => 0.1,						# initial sd of shock on hottest chain
		"past_iterations" => 30,						# num of periods used to compute Cov(p)
		"min_accept_tol"  => 100000,					# ABC-MCMC cutoff for rejecting small improvements
		"max_accept_tol"  => 100000,					# ABC-MCMC cutoff for rejecting small improvements
		"min_disttol"     => 0.1,						# distance tol for jumps from coldest chain
		"max_disttol"     => 0.1,						# distance tol for jumps from hottest chain
		"min_jump_prob"   => 0.05,						# prob of jumps from coldest chain
		"max_jump_prob"   => 0.05)						# prob of jumps from hottest chain
	# setup the BGP algorithm
	MA = MAlgoBGP(mprob,opts)
	# run the estimation
	runMOpt!(MA)
	fig = figure("parameter histograms") 
	plt[:hist](convert(Array,MOpt.parameters(MA.MChains[1])),15)
	return MA
end