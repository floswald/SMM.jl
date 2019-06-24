
module TestEconometrics

using FactCheck, DataFrames, MOpt, Lazy, DataFramesMeta

facts("checking standard errors") do 

	pb   = [ "p1" => [0.3, -1,1] , "p2" => [-0.3,-1,1]]
	moms = DataFrame(name=["mu1","mu2"],value=[0.0,0.0],weight=rand(2))
	mprob = @> MProb() addSampledParam!(pb) addMoment!(moms) addEvalFunc!(MOpt.objfunc_norm)
	opts =["N"=>3,"shock_var"=>1.0,"mode"=>"serial","maxiter"=>500,"path"=>".","maxtemp"=>100,"min_shock_sd"=>0.1,"max_shock_sd"=>1.0,"past_iterations"=>30,"min_disttol"=>0.1,"max_disttol"=>1.0,"min_jump_prob"=>0.05,"max_jump_prob"=>0.1] 

	MA = MAlgoBGP(mprob,opts)
	run!(MA)

	m = MOpt.score(MA)
	println(m)

	se = std(MA)
	println(se)

	println("se")

	# get 5% and 95% quantiles of param distributions
	chs = MOpt.parameters_ID(MA.MChains)
	ch1 = chs[chs[:chain_id].==1,convert(Array{Symbol}, sort(collect(keys(MA.m.params_to_sample))))]

	# the 5 and 95 quantile of the param distribution should correspond to
	# q_05 = mu - 1.96 * sigma/sqrt(n), where sigma is the estimated se
	# q_96 = mu + 1.96 * sigma/sqrt(n), where sigma is the estimated se

	mu = @> begin
       ch1
       @select(p1 = mean(:p1),p2=mean(:p2))
    end


    qtiles = Dict()
    for p in names(ch1)
    	qtiles[p] = quantile(ch1[p].data,[0.05,0.95])
    	println("empirical_05 = $(qtiles[p][1])")
    	println("MCMC_05 = $(mu[p][1] .- 1.96 * se[p]/sqrt(opts["N"]))")
    	@fact qtiles[p][1] => roughly(mu[p][1] .- 1.96 * se[p]/sqrt(opts["N"]),atol=0.05)
    	println("empirical_95 = $(qtiles[p][2])")
    	println("MCMC_95 = $(mu[p][1] .+ 1.96 * se[p]/sqrt(opts["N"]))")
    	@fact qtiles[p][2] => roughly(mu[p][1] .+ 1.96 * se[p]/sqrt(opts["N"]),atol=0.05)
    end

end


end
