
# The BGP MCMC Algorithm: Likelihood-Free Parallel Tempering
# ==========================================================
#
# http://link.springer.com/article/10.1007%2Fs11222-012-9328-6
# http://fr.arxiv.org/abs/1108.3423
#
# Baragatti, Grimaud and Pommeret (BGP)
# 
# Approximate Bayesian Computational (ABC) methods (or likelihood-free methods) have appeared in the past fifteen years as useful methods to perform Bayesian analyses when the likelihood is analytically or computationally intractable. Several ABC methods have been proposed: Monte Carlo Markov Chains (MCMC) methods have been developped by Marjoramet al. (2003) and by Bortotet al. (2007) for instance, and sequential methods have been proposed among others by Sissonet al. (2007), Beaumont et al. (2009) and Del Moral et al. (2009). Until now, while ABC-MCMC methods remain the reference, sequential ABC methods have appeared to outperforms them (see for example McKinley et al. (2009) or Sisson et al. (2007)). In this paper a new algorithm combining population-based MCMC methods with ABC requirements is proposed, using an analogy with the Parallel Tempering algorithm (Geyer, 1991). Performances are compared with existing ABC algorithms on simulations and on a real example.
export jumpParams!

# Define a Chain Type for BGP
type MHChain <: Chain
	id::Int             # chain id
	i::Int              # current index
	infos        ::DataFrame   # DataFrameionary of arrays(L,1) with eval, ACC and others
	parameters   ::DataFrame   # DataFrameionary of arrays(L,1), 1 for each parameter
	moments      ::DataFrame   # DataFrameionary of DataArrays(L,1), 1 for each moment
	accept_tol   ::Float64     # acceptance tolerance: new is not a "substantial" improvement over old, don't accept
	dist_tol     ::Float64 # percentage value distance from current chain i that is considered "close" enough. i.e. close = ((val_i - val_j)/val_j < dist_tol )

	params_nms   ::Array{Symbol,1}	# names of parameters (i.e. exclusive of "id" or "iter", etc)
	moments_nms  ::Array{Symbol,1}	# names of moments
	params2s_nms ::Array{Symbol,1}  # DataFrame names of parameters to sample 

	# TODO need either of those not both
	# the paper uses tempering to set up the kernel,
	# tibo uses shock_sd to amplify the shocks. expect small difference.
	shock_sd   ::Float64 # sd of shock to 

	function MHChain(id,MProb,L,temp,shock,accept_tol,dist_tol,jump_prob)
		this       = Chain(MProb,L)
		return this
    end
end

type MHChains
	MChains :: Array{MHChain,1}
end

type MAlgoMH <: MAlgo
    m               :: MProb # an MProb
    opts            :: Dict	# list of options
    i               :: Int 	# iteration
    current_param   :: Array{Dict,1}  # current param value: one Dict for each chain
    MChains         :: Array{BGPChain,1} 	# collection of Chains: if N==1, length(chains) = 1
  
    function MAlgoMH(m::MProb,opts=["N"=>3,"min_shock_sd"=>0.1,"max_shock_sd"=>1.0,"maxiter"=>100,"maxtemp"=> 100])

    	this      = new()
    	this.m    = m 
    	this.opts = opts
    	this.i    = 0
    	this.current_param = [ deepcopy(m.initial_value) for i=1:opts["N"] ] 

		acctol       = linspace(opts["min_accept_tol"],opts["max_accept_tol"],opts["N"])
		shocksd      = linspace(opts["min_shock_sd"],opts["max_shock_sd"],opts["N"])
		disttol      = linspace(opts["min_disttol"],opts["max_disttol"],opts["N"])
		jump_prob    = linspace(opts["min_jump_prob"],opts["max_jump_prob"],opts["N"])
	  	this.MChains = [BGPChain(i,m,opts["maxiter"],temps[i],shocksd[i],acctol[i],disttol[i],jump_prob[i]) for i=1:opts["N"] ]

	  	return this
    end
end

# ---------------------------  BGP ALGORITHM ----------------------------

# computes new candidate vectors for each chain
# accepts/rejects that vector on each chain, according to some rule
# *) computes N new parameter vectors
# *) applies a criterion to accept/reject any new params
# *) stores the result in chains
function computeNextIteration!( algo::MAlgoBGP )
    # here is the meat of your algorithm:
    # how to go from p(t) to p(t+1) ?


    # check if we reached end of chain
	if algo.i == algo["maxiter"]
	    println("reached end of chain. goodbye.")
	    return true
	else
		# else update iteration count on all chains
		updateIterChain!(algo.MChains)

		# check algo index is the same on all chains
		for ic in 1:algo["N"]
			@assert algo.i == algo.MChains[ic].i
		end

		# New Candidates
		# --------------
		if algo.i > 1
			# MVN = getParamKernel(algo)	# returns a MvNormal object
			MVN = getParamCovariance(algo)	# returns a Cov matrix
			getNewCandidates!(algo,MVN)
		end

		# evaluate objective on all chains
		# --------------------------------
		v = pmap( x -> evaluateObjective(algo.m,x), algo.current_param)


		# Part 1) LOCAL MOVES ABC-MCMC for i={1,...,N}. accept/reject
		# -----------------------------------------------------------
		localMovesMCMC!(algo,v)

		# Part 2) EXCHANGE MOVES 
		# ----------------------
		# starting mixing in period 4
		if algo.i>2 && algo["N"] > 1 
			exchangeMoves!(algo)
		end

		# Part 3) update sampling variances
	end
end

# notice: higher tempering draws candiates further spread out,
# but accepts lower function values with lower probability
function localMovesMCMC!(algo::MAlgoBGP,v::Array)
	for ch in 1:algo["N"]
		chain = algo.MChains[ch]
		eval_new  = v[ch]

		xold = -99.0
		if algo.i == 1
			# accept all
			ACC = true; prob = 1.0
			chain.infos[algo.i,:accept_rate] = 0.1
			eval_new.status = 1
			appendEval!(chain,eval_new,ACC,prob)
		else
			eval_old = getEval(chain,algo.i-1)

			prob = minimum([1.0, exp(chain.tempering *(eval_old.value - eval_new.value))])
			prob = prob * (eval_new.value < chain.accept_tol)

			if isna(prob)
				prob = 0.0
				eval_old.status = -1
				ACC = false
				appendEval!(chain,eval_old,ACC,prob)

			elseif !isfinite(xold)
				prob = 1.0
				eval_new.status = -2
				ACC = false
				appendEval!(chain,eval_new,ACC,prob)
			else 
				status = 1
				if prob > rand()
					ACC = true
					appendEval!(chain,eval_new,ACC,prob)
				else
					ACC = false
					appendEval!(chain,eval_old,ACC,prob)
				end
			end 
		    # update sampling variances
		    chain.infos[algo.i,:accept_rate]   = 0.9 * chain.infos[algo.i-1,:accept_rate] + 0.1 * ACC
		    chain.shock_sd                     = chain.shock_sd * (1+ 0.05*( 2*(chain.infos[algo.i,:accept_rate]>0.234) -1) )
		    chain.infos[algo.i,:shock_sd]      = chain.shock_sd
		    chain.infos[algo.i,:perc_new_old] = (eval_new.value - eval_old.value) / abs(eval_old.value)
		end
	end
end


# function getNewCandidates!(algo::MAlgoBGP,MVN::MvNormal)
function getNewCandidates!(algo::MAlgoBGP,VV::Matrix)

	# update chain by chain
	for ch in 1:algo["N"]

		# TODO
		# getParamKernel could be in here
		# chain.temperature should parameterize MVN somehow (as in their toy example: multiplies the variance)
		# setup a MvNormal
		# VV2 = VV.*algo.MChains[ch].tempering
		VV2 = VV
		MVN = MvNormal(VV2)

		# constraint the shock_sd: 95% conf interval should not exceed overall param interval width
		shock_ub  = minimum(  (algo.m.params_to_sample_df[:ub] - algo.m.params_to_sample_df[:lb]) ./ (1.96 * 2 *diag(VV2) ))
		algo.MChains[ch].shock_sd  = min(algo.MChains[ch].shock_sd , shock_ub)

		# shock parameters on chain index ch
		shock = rand(MVN) * algo.MChains[ch].shock_sd

		if get(algo.opts,"print_level",0) > 3
			Base.info("shock to parameters on chain $ch :")
			Base.info(shock)
		end

		updateCandidateParam!(algo,ch,shock)

	end

end


function jumpParams!(algo::MAlgoBGP,ch::Int,shock::Dict)
	eval_old = getLastEval(algo.MChains[ch])
	for k in keys(eval_old.params)
		algo.current_param[ch][k] = fitMirror(eval_old.params[k] + shock[k] , 
												algo.m.params_to_sample[k][:lb],
												algo.m.params_to_sample[k][:ub])
	end
end


# save algo chains component-wise to HDF5 file
function save(algo::MAlgoBGP, filename::ASCIIString)
    # step 1, create the file if it does not exist
    ff5 = h5open(filename, "w")

	# saving the opts dict: complicated because values are numbers and strings.
    vals = ASCIIString[]
    keys = ASCIIString[]
	for (k,v) in algo.opts
		if typeof(v) <: Number
			push!(vals,"$v")
		else
			push!(vals,v)
		end
		push!(keys,k)
	end
    write(ff5,"algo/opts/keys",keys)
    write(ff5,"algo/opts/vals",vals)

	# saving the chains
	for cc in 1:algo["N"]
	    saveChainToHDF5(algo.MChains[cc], ff5, "chain/$cc")
	end

    close(ff5)
end


function show(io::IO,MA::MAlgoBGP)
	print(io,"BGP Algorithm with $(MA["N"]) chains\n")
	print(io,"====================================\n")
	print(io,"\n")
	print(io,"Algorithm\n")
	print(io,"---------\n")
	print(io,"Current iteration: $(MA.i)\n")
	print(io,"Number of params to estimate: $(length(MA.m.p2sample_sym))\n")
	print(io,"Number of moments to match: $(length(MA.m.moments_subset))\n")
	print(io,"Chains\n")
	print(io,"------\n")
	print(io,"Tempering range: [$(MA.MChains[1].tempering),$(MA.MChains[end].tempering)]\n")
	print(io,"Jump probability range: [$(MA.MChains[1].jump_prob),$(MA.MChains[end].jump_prob)]\n")
end
