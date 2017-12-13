	

function parallelNormal(niter=200)
	# data are generated from a bivariate normal
	# with mu = [a,b] = [0,0]
	# aim: 
	# 1) sample [a',b'] from a space [-3,3] x [-2,2] and
	# 2) find true [a,b] by computing distance(S([a',b']), S([a,b]))
	#    and accepting/rejecting [a',b'] according to BGP
	# 3) S([a,b]) returns a summary of features of the data: 2 means

	# initial value
	pb    = OrderedDict("p1" => [0.2,-3,3] , "p2" => [-0.2,-2,2] )
	moms = DataFrame(name=["mu1","mu2"],value=[-1.0,1.0],weight=ones(2))
	mprob = MProb() 
	addSampledParam!(mprob,pb) 
	addMoment!(mprob,moms) 
	addEvalFunc!(mprob,objfunc_norm)

	nchains = 3

	opts =Dict("N"=>nchains,
		"maxiter"=>niter,
		"maxtemp"=> 5,
            # choose inital sd for each parameter p
            # such that Pr( x \in [init-b,init+b]) = 0.975
            # where b = (p[:ub]-p[:lb])*opts["coverage"] i.e. the fraction of the search interval you want to search around the initial value
		"coverage"=>0.025,  # i.e. this gives you a 95% CI about the current parameter on chain number 1.
		"sigma_update_steps"=>10,
		"sigma_adjust_by"=>0.01,
		"smpl_iters"=>1000,
		"parallel"=>true,
		"maxdists"=>[0.05 for i in 1:nchains],
		"mixprob"=>0.3,
		"acc_tuner"=>12.0,
		"animate"=>false)

	# plot slices of objective function
	s = doSlices(mprob,30)
	plot(s,:value)  # plot objective function over param values
	savefig(joinpath(dirname(@__FILE__),"../../slices-v.png"))
	plot(s,:mu1)  # plot value of moment :mu1 over param values
	savefig(joinpath(dirname(@__FILE__),"../../slices-m.png"))
	plot(s,:mu2)  # plot value of moment :mu2 over param values
	savefig(joinpath(dirname(@__FILE__),"../../slices-m2.png"))

	# setup the BGP algorithm
	MA = MAlgoBGP(mprob,opts)

	# run the estimation
	runMOpt!(MA)
	@show summary(MA)

	histogram(MA.chains[1]);
	savefig(joinpath(dirname(@__FILE__),"../../histogram.png"))
	plot(MA.chains[1]);
	savefig(joinpath(dirname(@__FILE__),"../../lines.png"))
	return MA
end

function serialNormal(niter=200)
	# data are generated from a bivariate normal
	# with mu = [a,b] = [0,0]
	# aim: 
	# 1) sample [a',b'] from a space [-3,3] x [-2,2] and
	# 2) find true [a,b] by computing distance(S([a',b']), S([a,b]))
	#    and accepting/rejecting [a',b'] according to BGP
	# 3) S([a,b]) returns a summary of features of the data: 2 means

	# initial value
	pb    = OrderedDict("p1" => [0.2,-3,3] , "p2" => [-0.2,-2,2] )
	moms = DataFrame(name=["mu1","mu2"],value=[-1.0,1.0],weight=ones(2))
	mprob = MProb() 
	addSampledParam!(mprob,pb) 
	addMoment!(mprob,moms) 
	addEvalFunc!(mprob,objfunc_norm)

	nchains = 3

	opts =Dict("N"=>nchains,
		"maxiter"=>niter,
		"maxtemp"=> 5,
            # choose inital sd for each parameter p
            # such that Pr( x \in [init-b,init+b]) = 0.975
            # where b = (p[:ub]-p[:lb])*opts["coverage"] i.e. the fraction of the search interval you want to search around the initial value
		"coverage"=>0.025,  # i.e. this gives you a 95% CI about the current parameter on chain number 1.
		"sigma_update_steps"=>10,
		"sigma_adjust_by"=>0.01,
		"smpl_iters"=>1000,
		"parallel"=>false,
		"maxdists"=>[0.05 for i in 1:nchains],
		"mixprob"=>0.3,
		"acc_tuner"=>12.0,
		"animate"=>false)

	# plot slices of objective function
	s = doSlices(mprob,30)
	plot(s,:value)  # plot objective function over param values
	savefig(joinpath(dirname(@__FILE__),"../../slices-v.png"))
	plot(s,:mu1)  # plot value of moment :mu1 over param values
	savefig(joinpath(dirname(@__FILE__),"../../slices-m.png"))
	plot(s,:mu2)  # plot value of moment :mu2 over param values
	savefig(joinpath(dirname(@__FILE__),"../../slices-m2.png"))

	# setup the BGP algorithm
	MA = MAlgoBGP(mprob,opts)

	# run the estimation
	runMOpt!(MA)
	@show summary(MA)

	histogram(MA.chains[1]);
	savefig(joinpath(dirname(@__FILE__),"../../histogram.png"))
	plot(MA.chains[1]);
	savefig(joinpath(dirname(@__FILE__),"../../lines.png"))
	return MA
end


function BGP_example()

	p = OrderedDict("theta" => [2.0,-2,10])
	mprob = MProb()
	addSampledParam!(mprob,p)
	addEvalFunc!(mprob,MomentOpt.objfunc_BGP)

	nchains = 15

	opts =Dict("N"=>nchains,
		"maxiter"=>500,
		"maxtemp"=> 4,
            # choose inital sd for each parameter p
            # such that Pr( x \in [init-b,init+b]) = 0.975
            # where b = (p[:ub]-p[:lb])*opts["coverage"] i.e. the fraction of the search interval you want to search around the initial value
		"coverage"=>0.005,  # i.e. this gives you a 95% CI about the current parameter on chain number 1.
		"maxdists"=>linspace(0.025, 2,nchains),
		"mixprob"=>0.3,
		"acc_tuner"=>12.0,
		"animate"=>false)
	# setup the BGP algorithm
	MA = MAlgoBGP(mprob,opts)
	# run the estimation
	runMOpt!(MA)
	@show summary(MA)

	histogram(MA.chains[1])
end

function serialSlices()
	# data are generated from a bivariate normal
	# with mu = [a,b] = [0,0]
	# aim: 
	# 1) sample [a',b'] from a space [-1,1] x [-1,1] and
	# 2) find true [a,b] by computing distance(S([a',b']), S([a,b]))
	#    and accepting/rejecting [a',b'] according to BGP
	# 3) S([a,b]) returns a summary of features of the data

	# initial value
	pb    = OrderedDict("p1" => [-1.0,-3,3] , "p2" => [1.0,-2,2] )
	moms = DataFrame(name=["mu1","mu2"],value=[-1.0,1.0],weight=ones(2))
	mprob = MProb() 
	addSampledParam!(mprob,pb) 
	addMoment!(mprob,moms) 
	addEvalFunc!(mprob,objfunc_norm)

	s = doSlices(mprob,50)
	p =plot(s,:value)
	display(p)
	return s
end
