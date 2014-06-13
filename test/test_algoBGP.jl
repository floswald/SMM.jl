


module TestAlgoBGP

using FactCheck, DataFrames

include("../src/mopt.jl")


p    = ["a" => 3.1 , "b" => 4.9]
pb   = [ "a" => [0,1] , "b" => [0,1] ]
moms = [
	"alpha" => [ 0.8 , 0.02 ],
	"beta"  => [ 0.8 , 0.02 ],
	"gamma" => [ 0.8 , 0.02 ]
]

mprob = Mopt.MProb(p,pb,Mopt.Testobj,moms)

facts("testing MAlgoBGP Constructor") do

	context("checking members") do

		opts =["N"=>5,"shock_var"=>1.0,"mode"=>"serial","maxiter"=>100,"path"=>".","maxtemp"=>100,"min_shock_sd"=>0.1,"max_shock_sd"=>1.0,"past_iterations"=>30,"min_jumptol"=>0.1,"max_jumptol"=>1.0] 

		MA = Mopt.MAlgoBGP(mprob,opts)

		@fact isa(MA,Mopt.MAlgo) => true
		@fact isa(MA,Mopt.MAlgoBGP) => true
		@fact isa(MA.m,Mopt.MProb) => true

		@fact MA.i => 0
		@fact length(MA.MChains) => opts["N"]

		for ix = 1:opts["N"]
			@fact MA.current_param[ix] => p 
		end

	end

	context("checking getters/setters on MAlgo") do

		opts =["N"=>5,"shock_var"=>1.0,"mode"=>"serial","maxiter"=>100,"path"=>".","maxtemp"=>100,"min_shock_sd"=>0.1,"max_shock_sd"=>1.0,"past_iterations"=>30,"min_jumptol"=>0.1,"max_jumptol"=>1.0] 
		MA = Mopt.MAlgoBGP(mprob,opts)

		# getters
		@fact MA["N"] => opts["N"]
		@fact MA["mode"] => opts["mode"]

		# setter
		MA["newOption"] = "hi"
		@fact MA["newOption"] => "hi"

	end

end



facts("testing getNewCandidates(MAlgoBGP)") do


	context("test checkbounds!(df,dict)") do

		df = DataFrame(a=8,b=-10)
		Mopt.checkbounds!(df,pb)
		@fact df[1,:a] => pb["a"][2]
		@fact df[1,:b] => pb["b"][1]

	end

	context("test getParamKernel") do

		opts =["N"=>5,"shock_var"=>1.0,"mode"=>"serial","maxiter"=>100,"path"=>".","maxtemp"=>100,"min_shock_sd"=>0.1,"max_shock_sd"=>1.0,"past_iterations"=>30,"min_jumptol"=>0.1,"max_jumptol"=>1.0] 
		MA = Mopt.MAlgoBGP(mprob,opts)

		# fill chains with random values
		for iter =1:(MA["maxiter"]-1)
			# update all chains to index 1
			Mopt.updateIter!(MA.MChains)

			# set a parameter vector on all chains using appendEval!
			myp = ["a" => rand() , "b" => rand()]
			mym = ["alpha" => rand(),"beta"  => rand(),"gamma" =>  rand()]
			ret = ["value" => 1.1, "params" => myp, "time" => 0, "status" => 1, "moments" => mym]
			for ich in 1:MA["N"]
				Mopt.appendEval!(MA.MChains[ich],ret,true,1)
			end
		end

		# get all parameters
		pars = Mopt.Allparameters(MA.MChains)

		# get the last MA["past_iterations"] iterations from each chain
		# get parameter_to_sample names as symbols 
		
		sub = pars[pars[:iter] .<= maximum([1,MA.MChains[1].i-MA["past_iterations"]]),MA.m.p2sample_sym]

		# get covariance matrix of those
		VV = cov(Mopt.array(sub)) + 0.0001 * Diagonal([1 for i=1:length(p)])

		# get kernel and check VV
		myVV = Mopt.getParamKernel(MA)

		@fact cov(myVV)[:] == VV[:] => true

	end

	context("test updateCandidateParam!") do

		# taking a specific MvNormal instance, 
		# I can know exactly how parameters should change when shocked

		# I'm interested in how algo.candidate_param changes on chain ch

		# setup an Algo
		opts =["N"=>5,"shock_var"=>1.0,"mode"=>"serial","maxiter"=>100,"path"=>".","maxtemp"=>100,"min_shock_sd"=>0.1,"max_shock_sd"=>1.0,"past_iterations"=>30,"min_jumptol"=>0.1,"max_jumptol"=>1.0] 
		MA = Mopt.MAlgoBGP(mprob,opts)

		# fill chains with random values up to iteration ix
		ix = 5
		for iter =1:ix
			# update all chains to index 1
			Mopt.updateIter!(MA.MChains)

			# set a parameter vector on all chains using appendEval!
			for ich in 1:MA["N"]
				myp = ["a" => rand() , "b" => rand()]
				mym = ["alpha" => rand(),"beta"  => rand(),"gamma" =>  rand()]
				ret = ["value" => 1.1, "params" => myp, "time" => 0, "status" => 1, "moments" => mym]
				Mopt.appendEval!(MA.MChains[ich],ret,true,1)
			end
		end

		# get a "kernel"
		# pos def matrix:
		myVV = rand(length(MA.m.params_to_sample),length(MA.m.params_to_sample))
		myVV = myVV * myVV'
		MVN = Mopt.MvNormal(myVV)

		# manually create next parameter guess on
		# each chain
		newp = Array(DataFrame,MA["N"])
		shock = {i=>Array{Float64,1} for i=1:MA["N"]}

		# on all chains, the parameter entry number ix+1 
		# must correspond to entry number ix plus a shock from the kernel
		for ich in 1:MA["N"]
			shock[ich] = rand(MVN) * MA.MChains[ich].shock_sd 

			oldp = Mopt.parameters(MA.MChains[ich],ix-1,true)	# get a dataframe of row ix-1
			newp[ich] = copy(oldp[MA.m.p2sample_sym])	# get just those you wish to sample

			# add shock
			for i in 1:ncol(newp[ich])
				newp[ich][1,i] += shock[ich][i]
			end			

			# adjust for bounds
			Mopt.checkbounds!(newp[ich],MA.m.params_to_sample)
		end


		# check on algo.candidate_param
		for ich in 1:MA["N"]
			# call library funciton
			Mopt.updateCandidateParam!(MA,ich,shock[ich])

			@fact collect(values(MA.candidate_param[ich])) == map(x->x[1],collect(values(Mopt.df2dict(newp[ich])))) => true

		end


	end

	context("testing swapRows!") do

		opts =["N"=>5,"shock_var"=>1.0,"mode"=>"serial","maxiter"=>100,"path"=>".","maxtemp"=>100,"min_shock_sd"=>0.1,"max_shock_sd"=>1.0,"past_iterations"=>30,"min_jumptol"=>0.1,"max_jumptol"=>1.0] 
		MA = Mopt.MAlgoBGP(mprob,opts)

		# fill chains with random values up to iteration ix
		ix = 5
		for iter =1:ix
			# update all chains to index 1
			Mopt.updateIter!(MA.MChains)

			# set a parameter vector on all chains using appendEval!
			for ich in 1:MA["N"]
				myp = ["a" => rand() , "b" => rand()]
				mym = ["alpha" => rand(),"beta"  => rand(),"gamma" =>  rand()]
				ret = ["value" => 1.1, "params" => myp, "time" => 0, "status" => 1, "moments" => mym]
				Mopt.appendEval!(MA.MChains[ich],ret,true,1)
			end
		end

		# get a pair (i,j) of chains
		pair = Mopt.sample(MA.Jump_register,1)[1]

		# get params and moms
		p1 = Mopt.parameters(MA.MChains[pair[1]],ix)
		p2 = Mopt.parameters(MA.MChains[pair[2]],ix)
		m1 = Mopt.moments(MA.MChains[pair[1]],ix)
		m2 = Mopt.moments(MA.MChains[pair[2]],ix)


		# exchange
		Mopt.swapRows!(MA,pair,ix)

		# check parameters
		@fact Mopt.parameters(MA.MChains[pair[1]],ix) == p2 => true
		@fact Mopt.parameters(MA.MChains[pair[2]],ix) == p1 => true
		@fact Mopt.moments(MA.MChains[pair[1]],ix) == m2 => true
		@fact Mopt.moments(MA.MChains[pair[2]],ix) == m1 => true

	end

end


facts("testing MAlgo methods") do
	
	context("testing evaluateObjective(algo)") do

		opts =["N"=>5,"shock_var"=>1.0,"mode"=>"serial","maxiter"=>100,"path"=>".","maxtemp"=>100,"min_shock_sd"=>0.1,"max_shock_sd"=>1.0,"past_iterations"=>30,"min_jumptol"=>0.1,"max_jumptol"=>1.0] 
		MA = Mopt.MAlgoBGP(mprob,opts)

		which_chain = 1
		x = Mopt.evaluateObjective(MA,which_chain)

		@fact haskey(x,"value") => true
		@fact x["params"] => p
		@fact haskey(x,"moments") => true

		# change p on MA:
		newp = ["a" => 103.1 , "b" => -2.2]
		MA.candidate_param[which_chain] = newp
		x = Mopt.evaluateObjective(MA,1)
		@fact x["params"] => newp
	end


	context("updateChains!(algo)") do

		opts =["N"=>5,"shock_var"=>1.0,"mode"=>"serial","maxiter"=>100,"path"=>".","maxtemp"=>100,"min_shock_sd"=>0.1,"max_shock_sd"=>1.0,"past_iterations"=>30,"min_jumptol"=>0.1,"max_jumptol"=>1.0] 
		MA = Mopt.MAlgoBGP(mprob,opts)

		@fact MA.MChains[1].i => 0
		@fact MA.i => 0

		Mopt.updateChains!(MA)

		# compute objective function once
		v = Mopt.Testobj(p,moms,["alpha", "beta","gamma"])

		# each chain must have those values
		for ix in 1:opts["N"]
			ch = Mopt.getchain(MA,ix)	# get each chain
			@fact ch.i => 1
			@fact Mopt.evals(ch,1)[1] => v["value"]
			@fact Mopt.accept(ch,1)[1] => true # all true in first eval
			@fact map(x->x[1],collect(values(Mopt.parameters(ch,1)))) == collect(values(v["params"])) => true
			@fact map(x->x[1],collect(values(Mopt.moments(ch,1)))) == collect(values(v["moments"])) => true
		end # all chains
	end
end





end # module

