


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

		opts =["N"=>3,"shock_var"=>1.0,"mode"=>"serial","maxiter"=>100,"path"=>".","maxtemp"=>100,"min_shock_sd"=>0.1,"max_shock_sd"=>1.0] 

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

		opts =["N"=>3,"shock_var"=>1.0,"mode"=>"serial","maxiter"=>100,"path"=>".","maxtemp"=>100,"min_shock_sd"=>0.1,"max_shock_sd"=>1.0] 
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

		opts =["N"=>3,"shock_var"=>1.0,"mode"=>"serial","maxiter"=>100,"path"=>".","maxtemp"=>100,"min_shock_sd"=>0.1,"max_shock_sd"=>1.0,"past_iterations"=>30] 
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
		par2sample_name   = collect(keys(MA.m.params_to_sample))
		par2sample_sym    = Array(Symbol,length(par2sample_name))
		for i in 1:length(par2sample_name)
			par2sample_sym[i] = symbol(par2sample_name[i])
		end
		sub = pars[pars[:iter] .<= maximum([1,MA.MChains[1].i-MA["past_iterations"]]),par2sample_sym]

		# get covariance matrix of those
		VV = cov(Mopt.array(sub)) + 0.0001 * Diagonal([1 for i=1:length(p)])

		# get kernel and check VV
		myVV = Mopt.getParamKernel(MA)

		@fact cov(myVV)[:] == VV[:] => true

	end
end


facts("testing MAlgo methods") do
	
	context("testing evaluateObjective(algo)") do

		opts =["N"=>3,"shock_var"=>1.0,"mode"=>"serial","maxiter"=>100,"path"=>".","maxtemp"=>100,"min_shock_sd"=>0.1,"max_shock_sd"=>1.0] 
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

		opts =["N"=>3,"shock_var"=>1.0,"mode"=>"serial","maxiter"=>100,"path"=>".","maxtemp"=>100,"min_shock_sd"=>0.1,"max_shock_sd"=>1.0] 
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

