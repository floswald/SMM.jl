
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
  jumptol    ::Float64 # tolerance of chain i for "closeness" to chain j
  tempering  ::Float64 # tempering in update probability
  shock_sd   ::Float64 # sd of shock to 

  function BGPChain(id,MProb,L,temp,shock,tol)
    infos      = { "evals" => @data([0.0 for i = 1:L]) , "accept" => @data([false for i = 1:L]), "status" => [0 for i = 1:L], "exchanged_with" => [0 for i = 1:L]}
    parameters = { x => zeros(L) for x in ps_names(MProb) }
    moments    = { x => @data([0.0 for i = 1:L]) for x in ms_names(MProb) }
    return new(id,0,infos,parameters,moments,tol,temp,shock)
  end
end



type MAlgoBGP <: MAlgo
  m               :: MProb # an MProb
  opts            :: Dict	# list of options
  i               :: Int 	# iteration
  current_param   :: Array{Dict,1}  # current param value: one Dict for each chain
  candidate_param :: Array{Dict,1}  # dict of candidate parameters: if rejected, go back to current
  MChains         :: Array{BGPChain,1} 	# collection of Chains: if N==1, length(chains) = 1
  Jump_register   :: Array{(Int,Int),1}

  function MAlgoBGP(m::MProb,opts=["N"=>3,"min_shock_sd"=>0.1,"max_shock_sd"=>1.0,"mode"=>"serial","maxiter"=>100,"maxtemp"=> 100])

  	# temperatures for each chain
	temps = linspace(1,opts["maxtemp"],opts["N"])
  	# shock standard deviations for each chain
	shocksd = linspace(opts["min_shock_sd"],opts["max_shock_sd"],opts["N"])
  	# standard deviations for each chain
	jumptol = linspace(opts["min_jumptol"],opts["max_jumptol"],opts["N"])
  	# create chains
  	chains = [BGPChain(i,m,opts["maxiter"],temps[i],shocksd[i],jumptol[i]) for i=1:opts["N"] ]
  	# current param values
  	cpar = [ m.initial_value for i=1:opts["N"] ] 
  	# jump register
  	Jreg = (Int,Int)[]
  	for i in 1:opts["N"]
	  	for j in i:opts["N"]
	  		if i!=j
	  			push!(Jreg,(i,j))
	  		end
	  	end
	end


    return new(m,opts,0,cpar,cpar,chains,Jreg)
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

function runMopt( algo::MAlgoBGP )

	for i in 1:algo["maxiter"]
		computeNextIteration!( algo )
	end
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
			MVN = getParamKernel(algo)	# returns a MvNormal object
			getNewCandidates!(algo,MVN)
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


		# Part 1) LOCAL MOVES ABC-MCMC for i={1,...,N}. accept/reject
		# ============================

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

		# Part 2a) EXCHANGE MOVES without rings
		# ======================
		# 1) N exchange moves are proposed
		# 2) N pairs are chosen uniformly among all possible pairs with replacment 
		# 3) exchange (z_i,theta_i) with (z_j,theta_j) if val_i - val_j < tol

		jump_pairs = sample(algo.Jump_register,algo["N"],replace=false)
		for pair in jump_pairs
			distance = abs(evals(algo.MChains[pair[1]],algo.MChains[pair[1]].i) - evals(algo.MChains[pair[2]],algo.MChains[pair[2]].i))
			if distance[1] < algo.MChains[pair[1]].jumptol
				swapRows!(algo,pair,algo.i)
			end
		end


		# Part 2b) EXCHANGE MOVES with rings
		# ======================
		# 1) N exchange moves are proposed
		# 2) a ring with at least two associated chains is chosen randomly from all rings
		#    a ring is a partition of tolerance space, i.e. R1 = [0,1], R2 = [1,1.5] etc
		# 3) exchange (z_i,theta_i) with (z_j,theta_j) if val_i - val_j < tol
		# for ch2 in 1:algo["N"]
		# 	if ch==ch2
		# 		#nothing
		# 	else

		# 	end
		# end

		# TODO need a function exchangeRowd!(algo,from,to)

	end
	@assert algo.i == algo.MChains[1].i
end

function swapRows!(algo::MAlgoBGP,pair::(Int,Int),i::Int)

	tmp1 = copy(alls(algo.MChains[pair[1]],i,true))
	# |-------|--------|--------|--------------|-------|----------|----------|-----------|----------|
	# | Row # | accept | status | exhanged_with| evals | beta     | a        | b         | alpha    |
	# | 1     | true   | 1      | 0            | 1.1   | 0.893083 | 0.163338 | 0.0192453 | 0.677312 |
	tmp1[:exchanged_with] = pair[1]

	tmp2 = copy(alls(algo.MChains[pair[2]],i,true))
	tmp2[:exchanged_with] = pair[2]

	# plug back into respective chains at correct index i
	fillinFields!(algo.MChains[pair[1]].parameters,tmp2,i)
	fillinFields!(algo.MChains[pair[1]].moments,tmp2,i)
	fillinFields!(algo.MChains[pair[1]].infos,tmp2,i)

	fillinFields!(algo.MChains[pair[2]].parameters,tmp1,i)
	fillinFields!(algo.MChains[pair[2]].moments,tmp1,i)
	fillinFields!(algo.MChains[pair[2]].infos,tmp1,i)
end


# get past parameter values to compute new candidates
# we draw new candidates for each chain from a joint normal
# that depends on parameter vectors on ALL chains
# pardf = Dataframe with stacked parameter df for each chain
function getParamKernel(algo::MAlgoBGP)


	pardf = Allparameters(algo.MChains)
	# |-------|----|------|----------|-----------|
	# | Row # | id | iter | a        | b         |
	# | 1     | 1  | 1    | 0.901685 | 0.983013  |
	# | 2     | 1  | 2    | 0.282522 | 0.763859  |
	# | 3     | 1  | 3    | 0.817773 | 0.268273  |
	# ...



	# select the last algo["past_iterations"] iterations from each chain
	lower_bound_index = maximum([1,algo.MChains[1].i-algo["past_iterations"]])


	# compute Var-Cov matrix of parameters_to_sample
	# plus some small random noise
	VV = cov(array(pardf[pardf[:iter].<=lower_bound_index, algo.m.p2sample_sym])) + 0.0001 * Diagonal([1 for i=1:length(algo.m.p2sample_sym)])

	# setup a MvNormal
	MVN = MvNormal(VV)
	return MVN
end


function getNewCandidates!(algo::MAlgoBGP,MVN::MvNormal)

	# update chain by chain
	for ch in 1:algo["N"]

		# TODO 
		# i think getParamKernel should be in here
		# chain.temperature should parameterize MVN somehow (as in their toy example: multiplies the variance)

		# shock parameters on chain index ch
		shock = rand(MVN) * algo.MChains[ch].shock_sd

		updateCandidateParam!(algo,ch,shock)

	end

end

function checkbounds!(df::DataFrame,di::Dict)
	if nrow(df) > 1
		error("can only process a single row")
	end
	dfbounds = collectFields(di,1:length(di),true)
	for c in names(df)
		if df[1,c] > dfbounds[2,c]
			df[1,c] = dfbounds[2,c]
		elseif df[1,c] < dfbounds[1,c]
			df[1,c] = dfbounds[1,c]
		end
	end
end


function updateCandidateParam!(algo::MAlgoBGP,ch::Int,shock::Array{Float64,1})

	# get last param on that particular chain
	# as dataframe row
	oldpar = parameters(algo.MChains[ch],algo.MChains[ch].i-1,true)

	# add shock to each column of newpar
	newpar = copy(oldpar[algo.m.p2sample_sym])
	for c in 1:ncol(newpar)
		newpar[1,c] += shock[c]
	end

	# do bounds checking on newpar
	checkbounds!(newpar,algo.m.params_to_sample)

	# set as dict on algo.current_param
	fillinFields!(algo.candidate_param[ch],newpar)
end



