	

function sliceOpt(tol)
	pb    = OrderedDict("p1" => [0.2,-3,3] , "p2" => [-0.2,-2,20] )
	moms = DataFrame(name=["mu1","mu2"],value=[-1.0,10.0],weight=ones(2))
	mprob = MProb() 
	addSampledParam!(mprob,pb) 
	addMoment!(mprob,moms) 
	addEvalFunc!(mprob,objfunc_norm)

	s = optSlices(mprob,30,tol=tol)
	return s
end

function mprob_ex()
	pb    = OrderedDict("p1" => [0.2,-3,3] , "p2" => [-0.2,-2,20] )
	moms = DataFrame(name=["mu1","mu2","var"],value=[-1.0,10.0,1.0],weight=rand(3))
	mprob = MProb() 
	addSampledParam!(mprob,pb) 
	addMoment!(mprob,moms) 
	addEvalFunc!(mprob,objfunc_norm2)
	return mprob
end


function score_moments()

	m = mprob_ex()
	MomentOpt.FD_gradient(m,Dict(m.initial_value))

end




function std_errors()

	# 0. given an mprob:
	m = mprob_ex()

	# 1. obtain a best parameter estimate p
	# s = optSlices(m,10,tol=5.0)
	# p = s[:best][:p]

	p = m.initial_value
	# 2. compute "Data" var-cov matrix Σ by generating H samples of simulated data using p 
	Σ = getSigma(m,p,300)

	# 3. compute score of moment function
	J = FD_gradient(m,p)

	# 4. put all together to get standard errors
	# first get weighting matrix 
	W = Diagonal([v[:weight] for (k,v) in m.moments])
	SE = pinv(J *  W * J') * (J * W * Σ * W * J') * pinv(J * W * J') 
	return Dict(zip(collect(keys(p)),diag(SE)))

end

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
		"min_improve"=>[0.05 for i in 1:nchains],
		"mixprob"=>0.3,
		"acc_tuners"=>[12.0 for i in 1:nchains],
		"animate"=>false)

	# plot slices of objective function
	# s = doSlices(mprob,30)
	# plot(s,:value)  # plot objective function over param values
	# savefig(joinpath(dirname(@__FILE__),"../../slices-v.png"))
	# plot(s,:mu1)  # plot value of moment :mu1 over param values
	# savefig(joinpath(dirname(@__FILE__),"../../slices-m.png"))
	# plot(s,:mu2)  # plot value of moment :mu2 over param values
	# savefig(joinpath(dirname(@__FILE__),"../../slices-m2.png"))

	# setup the BGP algorithm
	MA = MAlgoBGP(mprob,opts)

	# run the estimation
	runMOpt!(MA)
	@show summary(MA)

	# histogram(MA.chains[1]);
	# savefig(joinpath(dirname(@__FILE__),"../../histogram.png"))
	# plot(MA.chains[1]);
	# savefig(joinpath(dirname(@__FILE__),"../../lines.png"))
	return MA
end

function serialNormal(npars,niter;save=false)

	nchains = 3

	if npars == 2
		opts =Dict("N"=>nchains,
			"maxiter"=>niter,
			"maxtemp"=> 5,
	            # choose inital sd for each parameter p
	            # such that Pr( x \in [init-b,init+b]) = 0.975
	            # where b = (p[:ub]-p[:lb])*opts["coverage"] i.e. the fraction of the search interval you want to search around the initial value
			"coverage"=>0.02,   # how big should the step size of new params be?
			"smpl_iters"=>1000,
			"parallel"=>false,
			"min_improve"=>[0.0 for i in 1:nchains],
			"acc_tuners"=>[20;2;1.0],
			"animate"=>false)
	else
		opts =Dict("N"=>nchains,
		"maxiter"=>niter,
		"maxtemp"=> 15,
            # choose inital sd for each parameter p
            # such that Pr( x \in [init-b,init+b]) = 0.975
            # where b = (p[:ub]-p[:lb])*opts["coverage"] i.e. the fraction of the search interval you want to search around the initial value
		"coverage"=>0.05,   # how big should the step size of new params be?
		"smpl_iters"=>1000,
		"parallel"=>false,
		"min_improve"=>[0.0 for i in 1:nchains],
		"acc_tuners"=>[50.0;10;5],
		"animate"=>false)

	end
	o = snorm_impl(opts,niter,npar=npars,save=save)
	return o

end

function snorm_standard()
	# data are generated from a bivariate normal
	# with mu = [a,b] = [0,0]
	# aim: 
	# 1) sample [a',b'] from a space [-3,3] x [-2,2] and
	# 2) find true [a,b] by computing distance(S([a',b']), S([a,b]))
	#    and accepting/rejecting [a',b'] according to BGP
	# 3) S([a,b]) returns a summary of features of the data: 2 means

	# initial value

	pb = OrderedDict()
	pb["p1"] = [0.2,-3,3]
	pb["p2"] = [-0.2,-2,2]
	moms = DataFrame(name=["mu1","mu2"],value=[-1.0,1.0],weight=[1.0,1.0])

	nchains = 3

	opts =Dict("N"=>nchains,
		"maxiter"=>200,
		"maxtemp"=> 5,
            # choose inital sd for each parameter p
            # such that Pr( x \in [init-b,init+b]) = 0.975
            # where b = (p[:ub]-p[:lb])*opts["coverage"] i.e. the fraction of the search interval you want to search around the initial value
		"coverage"=>0.02,   # how big should the step size of new params be?
		"smpl_iters"=>1000,
		"parallel"=>false,
		"min_improve"=>[0.0 for i in 1:nchains],
		"acc_tuners"=>[20;2;1.0],
		"animate"=>false)
	info("These moments are our targets")
	info("Parameter p_i corresponds to moment m_i")
	println(moms)

	mprob = MProb() 
	addSampledParam!(mprob,pb) 
	MomentOpt.addMoment!(mprob,moms) 
	MomentOpt.addEvalFunc!(mprob,objfunc_norm)

	# setup the BGP algorithm
	MA = MAlgoBGP(mprob,opts)

	# run the estimation
	runMOpt!(MA)
	@show summary(MA)

	ph = histogram(MA.chains[1],nbins=19);
	savefig(joinpath(dirname(@__FILE__),"../../histogram0.png"))
	pp = plot(MA.chains[1]);
	savefig(joinpath(dirname(@__FILE__),"../../lines0.png")	)
	return (MA,ph,pp)
end



function snorm_6_taxi(tol;par=false)

	pb = OrderedDict()
	pb["p1"] = [0.2,-3,3]
	pb["p2"] = [-0.2,-2,2]
	pb["p3"] = [-0.3,-2,2]
	pb["p4"] = [-0.4,-2,2]
	pb["p5"] = [0.3,-2,2]
	pb["p6"] = [0.4,-2,2]
	moms = DataFrame(name=["mu1","mu2","mu3","mu4","mu5","mu6"],
		value=[-1.0,1.0,0.5,-0.5,0.7,-0.7],weight=ones(6))
	
	mprob = MProb() 
	addSampledParam!(mprob,pb) 
	MomentOpt.addMoment!(mprob,moms) 
	MomentOpt.addEvalFunc!(mprob,objfunc_norm)

	s = optSlices(mprob,3,tol=tol,parallel=par)
	return (moms,s)
	
end

function snorm_18(niter)
	# data are generated from a bivariate normal
	# with mu = [a,b] = [0,0]
	# aim: 
	# 1) sample [a',b'] from a space [-3,3] x [-2,2] and
	# 2) find true [a,b] by computing distance(S([a',b']), S([a,b]))
	#    and accepting/rejecting [a',b'] according to BGP
	# 3) S([a,b]) returns a summary of features of the data: 2 means

	# initial value

	pb = OrderedDict()
	pb["p1"] = [0.2,0,1]
	pb["p2"] = [0.002,0,0.02]
	pb["p3"] = [-0.3,-2,2]
	pb["p4"] = [-0.4,-2,2]
	pb["p5"] = [0.008,0,0.02]
	pb["p6"] = [0.4,-2,2]
	pb["p7"] = [0.2,0,1]
	pb["p8"] = [0.002,0,0.02]
	pb["p9"] = [-0.3,-2,2]
	pb["p10"] = [-0.4,-2,2]
	pb["p11"] = [0.008,0,0.02]
	pb["p12"] = [0.4,-2,2]
	pb["p13"] = [0.2,0,1]
	pb["p14"] = [0.002,0,0.02]
	pb["p15"] = [-0.3,-2,2]
	pb["p16"] = [-0.4,-2,2]
	pb["p17"] = [0.008,0,0.02]
	pb["p18"] = [0.4,-2,2]
	moms = DataFrame(name=["mu1",
						   "mu2",
						   "mu3",
						   "mu4",
						   "mu5",
						   "mu6",
						   "mu7",
						   "mu8",
						   "mu9",
						   "mu10",
						   "mu11",
						   "mu12",
						   "mu13",
						   "mu14",
						   "mu15",
						   "mu16",
						   "mu17",
						   "mu18"],
		value=[-1.0,0.0,0.5,-0.5,0.008,-0.7,-1.0,0.0,0.5,-0.5,0.008,-0.7,-1.0,0.0,0.5,-0.5,0.008,-0.7],weight=[1.0,0.1,0.5,0.5,0.008,0.7,1.0,0.1,0.5,0.5,0.008,0.7,1.0,0.1,0.5,0.5,0.008,0.7])

	nchains = 3

	opts =Dict("N"=>nchains,
		"maxiter"=>niter,
		"maxtemp"=> 5,
            # choose inital sd for each parameter p
            # such that Pr( x \in [init-b,init+b]) = 0.975
            # where b = (p[:ub]-p[:lb])*opts["coverage"] i.e. the fraction of the search interval you want to search around the initial value
		"coverage"=>0.02,   # how big should the step size of new params be?
		"smpl_iters"=>1000,
		"sigma"=>0.001,
		"parallel"=>false,
		"min_improve"=>[0.0 for i in 1:nchains],
		"acc_tuners"=>[20;2;1.0],
		"animate"=>false,
		"batch_size" => 1)
	info("These moments are our targets")
	info("Parameter p_i corresponds to moment m_i")
	println(moms)

	mprob = MProb() 
	addSampledParam!(mprob,pb) 
	MomentOpt.addMoment!(mprob,moms) 
	MomentOpt.addEvalFunc!(mprob,objfunc_norm)

	# setup the BGP algorithm
	MA = MAlgoBGP(mprob,opts)

	# run the estimation
	runMOpt!(MA)
	@show summary(MA)

	# ph = histogram(MA.chains[1],nbins=19);
	# savefig(joinpath(dirname(@__FILE__),"../../histogram6.png"))
	# pp = plot(MA.chains[1]);
	# savefig(joinpath(dirname(@__FILE__),"../../lines6.png")	)
	return MA
end

function snorm_standard6()
	# data are generated from a bivariate normal
	# with mu = [a,b] = [0,0]
	# aim: 
	# 1) sample [a',b'] from a space [-3,3] x [-2,2] and
	# 2) find true [a,b] by computing distance(S([a',b']), S([a,b]))
	#    and accepting/rejecting [a',b'] according to BGP
	# 3) S([a,b]) returns a summary of features of the data: 2 means

	# initial value

	pb = OrderedDict()
	pb["p1"] = [0.2,-3,3]
	pb["p2"] = [-0.2,-2,2]
	pb["p3"] = [-0.3,-2,2]
	pb["p4"] = [-0.4,-2,2]
	pb["p5"] = [0.3,-2,2]
	pb["p6"] = [0.4,-2,2]
	moms = DataFrame(name=["mu1","mu2","mu3","mu4","mu5","mu6"],
		value=[-1.0,1.0,0.5,-0.5,0.7,-0.7],weight=[1.0,1.0,0.5,0.5,0.7,0.7])

	nchains = 3

	opts =Dict("N"=>nchains,
		"maxiter"=>200,
		"maxtemp"=> 5,
            # choose inital sd for each parameter p
            # such that Pr( x \in [init-b,init+b]) = 0.975
            # where b = (p[:ub]-p[:lb])*opts["coverage"] i.e. the fraction of the search interval you want to search around the initial value
		"coverage"=>0.02,   # how big should the step size of new params be?
		"smpl_iters"=>1000,
		"parallel"=>false,
		"min_improve"=>[0.0 for i in 1:nchains],
		"acc_tuners"=>[20;2;1.0],
		"animate"=>false)
	info("These moments are our targets")
	info("Parameter p_i corresponds to moment m_i")
	println(moms)

	mprob = MProb() 
	addSampledParam!(mprob,pb) 
	MomentOpt.addMoment!(mprob,moms) 
	MomentOpt.addEvalFunc!(mprob,objfunc_norm)

	# setup the BGP algorithm
	MA = MAlgoBGP(mprob,opts)

	# run the estimation
	runMOpt!(MA)
	@show summary(MA)

	ph = histogram(MA.chains[1],nbins=19);
	savefig(joinpath(dirname(@__FILE__),"../../histogram6.png"))
	pp = plot(MA.chains[1]);
	savefig(joinpath(dirname(@__FILE__),"../../lines6.png")	)
	return (MA,ph,pp)
end
function snorm_impl(opts,niter=200;npar=2,save=false)
	# data are generated from a bivariate normal
	# with mu = [a,b] = [0,0]
	# aim: 
	# 1) sample [a',b'] from a space [-3,3] x [-2,2] and
	# 2) find true [a,b] by computing distance(S([a',b']), S([a,b]))
	#    and accepting/rejecting [a',b'] according to BGP
	# 3) S([a,b]) returns a summary of features of the data: 2 means

	# initial value

	mapRange(a1,a2,b1,b2,s) = b1 + (s-a1)*(b2-b1)/(a2-a1)

	pb = OrderedDict()
	pb["p1"] = [0.2,-3,3]
	pb["p2"] = [-0.2,-20,20]
	moms = DataFrame(name=["mu1","mu2"],value=[-1.0,10.0],weight=ones(2))
	srand(12)

	if npar > 2
		spaces = vcat(rand(2),4.0^(-4),reverse(linspace(0.25,0.45,npar-1).^(-4)))
		for j in 3:npar
			ii = j
			# spaces = rand()^(-4)
			pb["p$j"] = [mapRange(0,1,-spaces[ii],spaces[ii],rand()), -spaces[ii],spaces[ii]]
			# push!(moms, [Symbol("mu$j"); pb["p$j"][1]; pb["p$j"][1]])   # if moments are equal starting value - easy.
			# push!(moms, [Symbol("mu$j"); pb["p$j"][1] * 1.1; pb["p$j"][1]])   # if moments are 110% of starting value - easy.
			y = mapRange(0,1,-spaces[ii],spaces[ii],rand())
			push!(moms, [Symbol("mu$j"); y ; y])
		end
	end

	info("These moments are our targets")
	info("Parameter p_i corresponds to moment m_i")
	println(moms)

	mprob = MProb() 
	addSampledParam!(mprob,pb) 
	MomentOpt.addMoment!(mprob,moms) 
	MomentOpt.addEvalFunc!(mprob,objfunc_norm)



	# plot slices of objective function
	# s = doSlices(mprob,30)
	# plot(s,:value)  # plot objective function over param values
	# savefig(joinpath(dirname(@__FILE__),"../../slices-v.png"))
	# plot(s,:mu1)  # plot value of moment :mu1 over param values
	# savefig(joinpath(dirname(@__FILE__),"../../slices-m.png"))
	# plot(s,:mu2)  # plot value of moment :mu2 over param values
	# savefig(joinpath(dirname(@__FILE__),"../../slices-m2.png"))

	# setup the BGP algorithm
	MA = MAlgoBGP(mprob,opts)

	# run the estimation
	runMOpt!(MA)
	@show summary(MA)

	ph = histogram(MA.chains[1],nbins=19);
	if save
		savefig(joinpath(dirname(@__FILE__),"../../histogram.png"))
	end
	pp = plot(MA.chains[1]);
	if save
		savefig(joinpath(dirname(@__FILE__),"../../lines.png")	)
	end
	return (MA,ph,pp)
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
		"min_improve"=>linspace(0.025, 2,nchains),
		"mixprob"=>0.3,
		"acc_tuners"=>[12.0 for i in 1:nchains],
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
