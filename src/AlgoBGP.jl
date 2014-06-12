
# The BGP MCMC Algorithm: Likelihood-Free Parallel Tempering
# ==========================================================
#
# http://link.springer.com/article/10.1007%2Fs11222-012-9328-6
# http://fr.arxiv.org/abs/1108.3423
#
# Baragatti, Grimaud and Pommeret (BGP)
# 
# Approximate Bayesian Computational (ABC) methods (or likelihood-free methods) have appeared in the past fifteen years as useful methods to perform Bayesian analyses when the likelihood is analytically or computationally intractable. Several ABC methods have been proposed: Monte Carlo Markov Chains (MCMC) methods have been developped by Marjoramet al. (2003) and by Bortotet al. (2007) for instance, and sequential methods have been proposed among others by Sissonet al. (2007), Beaumont et al. (2009) and Del Moral et al. (2009). Until now, while ABC-MCMC methods remain the reference, sequential ABC methods have appeared to outperforms them (see for example McKinley et al. (2009) or Sisson et al. (2007)). In this paper a new algorithm combining population-based MCMC methods with ABC requirements is proposed, using an analogy with the Parallel Tempering algorithm (Geyer, 1991). Performances are compared with existing ABC algorithms on simulations and on a real example.


# Define a Chain Type for BGP
type BGPChain <: AbstractChain
  id::Int             # chain id
  i::Int              # current index
  infos      ::Dict   # dictionary of arrays(L,1) with eval, ACC and others
  parameters ::Dict   # dictionary of arrays(L,1), 1 for each parameter
  moments    ::Dict   # dictionary of DataArrays(L,1), 1 for each moment
  tempering  ::Float64
  shock_sd   ::Float64

  function BGPChain(id,MProb,L,temp,tol)
    infos      = { "evals" => @data([0.0 for i = 1:L]) , "accept" => @data([false for i = 1:L]), "status" => [0 for i = 1:L], }
    parameters = { x => zeros(L) for x in ps_names(MProb) }
    moments    = { x => @data([0.0 for i = 1:L]) for x in ms_names(MProb) }
    return new(id,0,infos,parameters,moments,temp,tol)
  end
end



type MAlgoBGP <: MAlgo
  m               :: MProb # an MProb
  opts            :: Dict	# list of options
  i               :: Int 	# iteration
  current_param   :: Array{Dict,1}  # current param value: one Dict for each chain
  candidate_param :: Array{Dict,1}  # dict of candidate parameters: if rejected, go back to current
  MChains         :: Array{BGPChain,1} 	# collection of Chains: if N==1, length(chains) = 1
  MVNormShock     :: MvNormal

  function MAlgoBGP(m::MProb,opts=["N"=>3,"min_shock_sd"=>0.1,"max_shock_sd"=>1.0,"mode"=>"serial","maxiter"=>100,"maxtemp"=> 100])

  	# temperatures for each chain
	temps = linspace(1,opts["maxtemp"],opts["N"])
  	# shock standard deviations for each chain
	shocksd = linspace(opts["min_shock_sd"],opts["max_shock_sd"],opts["N"])
  	# create chains
  	chains = [BGPChain(i,m,opts["maxiter"],temps[i],shocksd[i]) for i=1:opts["N"] ]
  	# current param values
  	cpar = [ m.initial_value for i=1:opts["N"] ] 

    return new(m,opts,0,cpar,cpar,chains)
  end
end

# function shockallp(p,shocksd,VV) 
#   # sh = rmultnorm(1,rep(0,nrow(VV)),VV) * shocksd
#   MN = MvNormal(VV)
#   sh = rand(MN) .* shocksd

#   # update value for each param
#   for (pp in colnames(sh)) {
#     p[[pp]] = fitMirror( p[[pp]] + sh[,pp] ,
#                               LB = cf$pdesc[pp,'lb'],
#                               UB = cf$pdesc[pp,'ub'])  
#   }

#   return p 
# end


function getchain( algo::MAlgoBGP, i::Int)
	algo.MChains[i]
end

# computes new candidate vectors for each chain
# accepts/rejects that vector on each chain, according to some rule
# *) computes N new parameter vectors
# *) applies a criterion to accept/reject any new params
# *) stores the result in chains
function computeNextIteration!( algo::MAlgoBGP  )
    # here is the meat of your algorithm:
    # how to go from p(t) to p(t+1) ?

    # check if we reached end of chain
	if algo.i == algo["maxiter"]
	    println("reached end of chain. goodbye.")
	    return true
	else
		# else update iteration count on all chains
		updateIter!(algo.MChains)


		# New Candidates
		# --------------

		if algo.i > 1
			getNewCandidates!(algo)
		end

		# evaluate objective on all chains
		# --------------------------------
		if algo["mode"] == "serial"
			v = map( x -> evaluateObjective(algo,x), 1:algo["N"])
		else
			v = pmap( x -> evaluateObjective(algo,x), 1:algo["N"])
		end

		# notice that at this point, v[chain_id]["param"] is the candidate param vector
		# we append v to the chain
		# so if candidate is rejected, go back into v[chain_id]["param"] and reset to previous param

		# accept/reject
		# =============

		for ch in 1:algo["N"]
			if algo.i == 1
				# accept all
				ACC = true
		  		algo.current_param[ch] = algo.candidate_param[ch] 
		  		status = 1
			else
				xold = evals(algo.MChains[ch],algo.i-1)
				xnew = v[ch]["value"]
				prob = minimum([1, exp(xold - xnew)])
				if isna(prob)
					prob = 0
					status = -1
				elseif !isfinite(xold)
					prob = 1
					status = -2
				else 
					if prob > rand()
						ACC = true
				  		algo.current_param[ch] = algo.candidate_param[ch] 
					else
						ACC = false
						v[ch]["params"] = algo.current_param[ch]	# reset param in output of obj to previous value
						v[ch]["moments"] = moments(algo.MChains[ch],algo.i-1)	# reset moments in output of obj to previous value
					end
					status = 1
				end
			end
		    # append values to MChains at index ch
		    appendEval!(algo.MChains[ch],v[ch],ACC,status)
		end
	end
end


function getNewCandidates!(algo::MAlgoBGP)

	# get past parameter values to compute new candidates
	# we draw new candidates for each chain from a joint normal
	# that depends on parameter vectors on ALL chains
	# pardf = Dataframe with stacked parameter df for each chain
	pardf = parameters(algo.MChains)
	# select the last 30 observations from each chain
	lower_bound_index = maximum(1,nrow(pardf)-30*length(algo.MChains))
	par2sample_sym     = Array{Symbol,1}
	par2sample_name = collect(keys(algo.m.params_to_sample))
	for i in 1:length(par2sample_name)
		par2sample_sym[i] = symbol(par2sample_name[i])
	end

	# compute Var-Cov matrix of parameters
	# plus some small random noise
	VV = cov(pardf[lower_bound_index, par2sample_sym]) + 0.0001 .* Diagonal([1 for i=1:length(par2sample_sym)])

	# setup a MvNormal
	MVN = MvNormal(VV)

	# update chain by chain
	for ch in 1:algo["N"]

		# get last param on that chain
		# as dataframe row
		oldpar = parameters(algo.MChains[ch],algo.i-1,true)

		# shock parameters on chain ch
		shock = rand(MVN) * algo.MChains[ch].shock_sd

		newpar = oldpar[par2sample_sym] .+ shock

		# set as dict on algo.current_param
		fillinFields!(algo.MChains[ch].candidate_param,newpar)

	end

end





