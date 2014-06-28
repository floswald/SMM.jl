

if banan 

	# minimize the rosenbrock function

	# global min is (1,1)
	p    = ["a" => 0.8 , "b" => 1.1]
	pb   = [ "a" => [-1.0,2.1] , "b" => [-1.0,2.1] ]
	moms = DataFrame(moment=["alpha"],data_value=[1.0],data_sd=rand(1))



	mprob = MProb(p,pb,MOpt.banana,moms)
	opts =[
		"N"=>6,
		"mode"=>"serial",
		"maxiter"=> 500,
		"path"=>".",
		"maxtemp"=>100,
		"min_shock_sd"=>0.1,
		"max_shock_sd"=>1,
		"past_iterations"=>30,
		"min_accept_tol"=>1000,
		"max_accept_tol"=>1000,
		"min_disttol"=>10.0,
		"max_disttol"=>10.0,
		"min_jump_prob"=>0.5,
		"max_jump_prob"=>0.5] 
	MA = MAlgoBGP(mprob,opts)

	runMopt!(MA)
	MOpt.figure(1)
	plot(MA,"acc")
	MOpt.figure(2)
	plot(MA,"params_time")
	MOpt.figure(3)
	plot(MA,"params_dist")

	
else




	# bivariate normal
	
	p    = ["a" => 0.9 , "b" => -0.9]
	pb   = [ "a" => [-1,1] , "b" => [-1,1] ]
	moms = DataFrame(moment=["alpha","beta"],data_value=[0.0,0.0],data_sd=rand(2))

	mprob = MProb(p,pb,MOpt.objfunc_norm2,moms)

	opts =[
		"N"=>6,
		"mode"=>"serial",
		"savefile" => joinpath(pwd(),"chains.jdl"),
		"maxiter"=> 500,
		"path"=>".",
		"maxtemp"=>100,
		"min_shock_sd"=>0.1,
		"max_shock_sd"=>1,
		"past_iterations"=>30,
		"min_accept_tol"=>100,
		"max_accept_tol"=>100,
		"min_disttol"=>0.1,
		"max_disttol"=>0.1,
		"min_jump_prob"=>0.05,
		"max_jump_prob"=>0.2] 

	MA = MAlgoBGP(mprob,opts)

	runMopt!(MA)
	MOpt.figure(1)
	plot(MA,"acc")
	MOpt.figure(2)
	plot(MA,"params_time")
	MOpt.figure(3)
	plot(MA,"params_dist")

	# save
	save(MA.MChains,MA["savefile"])
	

end
