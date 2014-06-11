


# here is a list of potential options we may want to set on each algo:

#   shock_var        :: Float64 # variance of shock
#   np_shock         :: Float64
#   save_freq        :: Int
#   N                :: Int   # number of chains
#   n_untempered     :: Int   # number of chains untemprered
#   maxiter          :: Int   # maximum number of iterations
#   run              :: Int   # number of run

#   # cluster setup
#   mode             :: ASCIIString # {'serial','mpi'}

#   # paths for I/O
#   path :: ASCIIString # "where/to/save/this"

#   prepared      :: Bool
#   current_param :: Dict   # current parameter value
#   chains        :: MChain   # object of type MChain
# end

type MAlgoRandom <: MAlgo
  m               :: MProb # an MProb
  opts            :: Dict	# list of options
  i               :: Int 	# iteration
  current_param   :: Array{Dict,1}  # current param value: one Dict for each chain
  candidate_param :: Array{Dict,1}  # dict of candidate parameters: if rejected, go back to current
  MChains         :: MChain 	# collection of Chains: if N==1, length(chains) = 1

  function MAlgoRandom(m::MProb,opts=["N"=>3,"shock_var"=>1.0,"mode"=>"serial","maxiter"=>100])

  	# create chains
  	chains = MChain(opts["N"],m,opts["maxiter"])
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


function getchain( algo::MAlgoRandom, i::Int)
	algo.MChains[i]
end

# computes new candidate vectors for each chain
# accepts/rejects that vector on each chain, according to some rule
# *) computes N new parameter vectors
# *) applies a criterion to accept/reject any new params
# *) stores the result in chains
function computeNextIteration!( algo::MAlgoRandom  )
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

		# # compute Var-Cov matrix
	 #    VV = cov(chains[lower_bound_index:nrow(chains),params_to_sample2]) + 0.0001 * diag(length(params_to_sample2))

		if algo.i > 1
			for ch in 1:algo["N"]
			  	# this updating rule can differ by chain!
			  	for p in algo.m.params_to_sample
			  		algo.candidate_param[ch][p] = algo.current_param[ch][p] + randn()*shock_var
			    end
			end
		else
			# candidate = initial_value, so ok for first iteration
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
		    appendEval!(algo.MChains,ch,v[ch],ACC,status)
		end
	end
end







