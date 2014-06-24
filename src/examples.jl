


if banan 

	# minimize the rosenbrock function

	# global min is (1,1)
	p    = ["a" => 0.8 , "b" => 1.1]
	pb   = [ "a" => [-1.0,2.1] , "b" => [-1.0,2.1] ]
	moms = [
		"m1" => [ 1.0 , 0.02 ],
		"m2"  => [ 1.0 , 0.02 ]
	]



	mprob = MOpt.MProb(p,pb,banana,moms)
	opts =["N"=>6,"mode"=>"serial","maxiter"=> 500,"path"=>".","maxtemp"=>100,"min_shock_sd"=>1,"max_shock_sd"=>3,"past_iterations"=>30,"min_accept_tol"=>500,"max_accept_tol"=>500,"min_disttol"=>100,"max_disttol"=>200,"min_jump_prob"=>0.05,"max_jump_prob"=>0.2] 
	MA = MOpt.MAlgoBGP(mprob,opts)

	MOpt.runMOpt(MA)

	x=MOpt.alls(MA.MChains)
	figure(1)
	for i in 1:6
		subplot(2,3,i)
		plot(x[x[:chain_id].==i,:a])
	end
	suptitle("chain estimates for a")
	figure(2)
	for i in 1:6
		subplot(2,3,i)
		plot(x[x[:chain_id].==i,:b])
	end
	suptitle("chain estimates for b")
else




	# bivariate normal
	
	p    = ["a" => 0.9 , "b" => -0.9]
	pb   = [ "a" => [-1,1] , "b" => [-1,1] ]
	moms = MOpt.DataFrame(moment=["alpha","beta"],data_value=[0.0,0.0],data_sd=rand(2))

	mprob = MOpt.MProb(p,pb,MOpt.objfunc_norm2,moms)

	opts =[
		"N"=>6,
		"mode"=>"serial",
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

	MA = MOpt.MAlgoBGP(mprob,opts)

	MOpt.runMopt!(MA)
	MOpt.figure(1)
	MOpt.plot(MA,"acc")
	MOpt.figure(2)
	MOpt.plot(MA,"params_time")
	MOpt.figure(3)
	MOpt.plot(MA,"params_dist")
	

end
