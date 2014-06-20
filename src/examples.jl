



if banan 

	# minimize the rosenbrock function
	# y = 100 * (b - a * a)^2 + (1 - a)^2
	function banana(x::Dict,mom::Dict,whichmom::Array{ASCIIString,1})

		model = 100 .* (x["b"] - x["a"].^2 ).^2 .+ (1.-x["a"])^2
		data  = 0.0

		value = mean((data .- model).^2)

		# we just want to find the lowest value - no moments involved.
		ret = ["value" => value, "params" => x, "time" => 1.0, "status" => 1, "moments" => mom]
		return ret
	end


	# global min is (1,1)
	p    = ["a" => 0.8 , "b" => 1.1]
	pb   = [ "a" => [-1.0,2.1] , "b" => [-1.0,2.1] ]
	moms = [
		"m1" => [ 1.0 , 0.02 ],
		"m2"  => [ 1.0 , 0.02 ]
	]



	mprob = Mopt.MProb(p,pb,banana,moms)
	opts =["N"=>6,"mode"=>"serial","maxiter"=> 500,"path"=>".","maxtemp"=>100,"min_shock_sd"=>1,"max_shock_sd"=>3,"past_iterations"=>30,"min_accept_tol"=>500,"max_accept_tol"=>500,"min_disttol"=>100,"max_disttol"=>200,"min_jump_prob"=>0.05,"max_jump_prob"=>0.2] 
	MA = Mopt.MAlgoBGP(mprob,opts)

	Mopt.runMopt(MA)


	x=Mopt.alls(MA.MChains)
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
	function objfunc_norm2(p::Dict,mom::Dict,whichmom::Array{ASCIIString,1})

		sigma = reshape([1.0, 0.0,0.0,1.0],2,2)
		ns = 5000

		# compute simulated moments
		mu = [p["a"],p["b"]]
		MVN = Mopt.MvNormal(mu,sigma) 

		# get data mometns
		muD = Float64[mom["m1"][1], mom["m2"][1]]

	    # simulate model moments 
	    moments = mean(rand(MVN,ns),2)

	    # value = data - model
		value = mean((muD - moments).^2)

		momdict = ["m1" => moments[1], "m2" => moments[2]]

		ret = ["value"=>value, "params" =>p, "time" =>0.0, "status" => 1, "moments" => momdict]
		return ret
	end

	p    = ["a" => 0.5 , "b" => -0.5]
	pb   = [ "a" => [-1,1] , "b" => [-1,1] ]
	moms = [
		"m1" => [ 0.0 , 0.02 ],
		"m2"  => [ 0.0 , 0.02 ]
	]

	mprob = Mopt.MProb(p,pb,objfunc_norm2,moms)

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

	MA = Mopt.MAlgoBGP(mprob,opts)

	Mopt.runMopt!(MA)
	Mopt.plot(MA,"acc")
	Mopt.plot(MA,"params_time")
	Mopt.plot(MA,"params_dist")
	

end
