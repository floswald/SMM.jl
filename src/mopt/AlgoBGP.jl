
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

type Chain
    evals     :: Array{Eval}
    id        :: Int64
    iter      :: Int64
    accepted  :: Array{Bool}
    exchanged :: Array{Int}
    m         :: Mprob
    sigmas    :: Vector{Float64}

    function Chain(id::Int,n::Int,m::Mprob,sig::Vector{Float64})
        this           = new()
        this.evals     = Array{Eval}(n)
        this.accepted  = falses(n)
        this.exchanged = zeros(Int,n)
        this.id        = id
        this.iter      = 1
        this.m         = m
        this.sigma     = PDiagMat(sig)
        return this
    end
end

function next_eval!(c::Chain)
    # generate new parameter vector from last accepted param
    pp = getNewCandidates(c)

    # evaluate objective 
    ev = evaluateObjective(c.m,pp)

    # accept reject 
    doAcceptRecject!(c,ev)


end

type MAlgoBGP <: MAlgo
    m               :: MProb # an MProb
    opts            :: Dict	# list of options
    i               :: Int 	# iteration
    current_param   :: Array{OrderedDict}  # current param value: one Dict for each chain
    MChains         :: Array{Chain} 	# collection of Chains: if N==1, length(chains) = 1
  
    function MAlgoBGP(m::MProb,opts=Dict("N"=>3,"maxiter"=>100,"maxtemp"=> 2,"bound_prob"=>0.15,"disttol"=>0.1))

        if opts["N"] > 1
    		temps     = linspace(1.0,opts["maxtemp"],opts["N"])
            # initial std dev for each parameter to achieve at least bound_prob on the bounds
            init_sd = OrderedDict( k => initvar(v[:lb],opts["bound_prob"]) for (k,v) in m.params_to_sample)
            chains = Chain[Chain(i,opts["maxiter"],collect(values(init_sd)) .* temps[i]) for i in 1:opts["N"]]
          else
            temps     = [1.0]
            init_sd = OrderedDict( k => initvar(v[:lb],opts["bound_prob"]) for (k,v) in m.params_to_sample)
            chains = Chain[Chain(i,opts["maxiter"],collect(values(init_sd)) .* temps[i]) for i in 1:opts["N"]]
        end
	  	# current param values
	  	cpar = OrderedDict[ deepcopy(m.initial_value) for i=1:opts["N"] ] 
	    return new(m,opts,0,cpar,chains)
    end
end

# this appends ACCEPTED values only.
function appendEval!(chain::BGPChain, ev:: Eval, ACC::Bool, prob::Float64)
    chain.infos[chain.i,:evals]  = ev.value
    chain.infos[chain.i,:prob]   = prob
    chain.infos[chain.i,:accept] = ACC
    chain.infos[chain.i,:status] = ev.status
    chain.infos[chain.i,:eval_time] = ev.time
    chain.infos[chain.i,:tempering] = chain.tempering
    for k in names(chain.moments)
        if !(k in [:chain_id,:iter])
            chain.moments[chain.i,k] = ev.simMoments[k]
        end
    end
    for k in names(chain.parameters)
        if !(k in [:chain_id,:iter])
            chain.parameters[chain.i,k] = ev.params[k]
        end
    end
    return nothing
end

# changing only the eval fields. used in swapRows!
function appendEval!(chain::BGPChain, ev:: Eval)
    chain.infos[chain.i,:evals]  = ev.value
    chain.infos[chain.i,:status] = ev.status
    chain.infos[chain.i,:eval_time] = ev.time
    for k in names(chain.moments)
        if !(k in [:chain_id,:iter])
            chain.moments[chain.i,k] = ev.simMoments[k]
        end
    end
    for k in names(chain.parameters)
        if !(k in [:chain_id,:iter])
            chain.parameters[chain.i,k] = ev.params[k]
        end
    end
    return nothing
end

# getEval(c::BGPChain,i::Int64) = getEval

# ---------------------------  BGP ALGORITHM ----------------------------

# computes new candidate vectors for each chain
# accepts/rejects that vector on each chain, according to some rule
# *) computes N new parameter vectors
# *) applies a criterion to accept/reject any new params
# *) stores the result in chains
function computeNextIteration!( algo::MAlgoBGP )
    # here is the meat of your algorithm:
    # how to go from p(t) to p(t+1) ?

	incrementChainIter!(algo.MChains)

	# check algo index is the same on all chains
	for ic in 1:algo["N"]
		@assert algo.i == algo.MChains[ic].i
	end

    # TODO 
    # this on each chain
    # START=========================================================
    pmap( x->next_eval!(x), algo.MChains ) # this does getNewCandidates, evaluateObjective, doAcceptRecject

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
	doAcceptRecject!(algo,v)

    # STOP=========================================================

	# Part 2) EXCHANGE MOVES only on master
	# ----------------------
	# starting mixing in period 3
	if algo.i>=2 && algo["N"] > 1 
		exchangeMoves!(algo)
	end
end

# notice: higher tempering draws candiates further spread out,
# but accepts lower function values with lower probability
function doAcceptRecject!(c::Chain,ev::Eval)
	if c.iter == 1
		# accept everything.
		ACC = true; prob = 1.0
		chain.infos[algo.i,:accept_rate] = 0.1
		eval_new.status = 1
		appendEval!(chain,eval_new,ACC,prob)
	else
		eval_old = getEval(chain,algo.i-1)

        if eval_new.status < 0
            prob = 0.0
            eval_old.status = -1
            ACC = false
            appendEval!(chain,eval_old,ACC,prob)
        else

            # prob = minimum([1.0,exp( chain.tempering *( eval_old.value - eval_new.value))])
			prob = minimum([1.0,exp( ( eval_old.value - eval_new.value) )])

			if isna(prob)
				prob = 0.0
				eval_old.status = -1
				ACC = false
				appendEval!(chain,eval_old,ACC,prob)

			elseif !isfinite(eval_old.value)
				prob = 1.0
				eval_new.status = -2
				ACC = false
				appendEval!(chain,eval_new,ACC,prob)
			else 
				# status = 1
				if prob > rand()
					ACC = true
					appendEval!(chain,eval_new,ACC,prob)
				else
					ACC = false
					appendEval!(chain,eval_old,ACC,prob)
				end
			end 
        end
	    # update sampling variances
        # -------------------------

        # update average acceptance rate: moving average
	    chain.infos[algo.i,:accept_rate]   = 0.9 * chain.infos[algo.i-1,:accept_rate] + 0.1 * ACC

        # update shock variance. want to achieve a long run accpetance rate of 23.4% (See Casella and Berger)
        accept_too_high = chain.infos[algo.i,:accept_rate]>0.234
        if accept_too_high
            chain.shock_sd *= 1.20  # increase variance by 10% => will accept less
        else # too low
            chain.shock_sd *= 0.9  # decrease variance by 10% => will accept more
        end
	    chain.infos[algo.i,:shock_sd]      = chain.shock_sd
	    chain.infos[algo.i,:perc_new_old] = (eval_new.value - eval_old.value) / abs(eval_old.value)
        debug("")
        debug("chain: $ch")
	    debug("ACCEPTED: $ACC")
        debug("accept_rate = $(chain.infos[algo.i,:accept_rate])")
	    debug("old value: $(eval_old.value)")
	    debug("new value: $(eval_new.value)")
        debug("")
	end
	end
end



function exchangeMoves!(algo::MAlgoBGP)
	# for all chains
	for ch in 1:algo["N"]
		e1 = getEval(algo.MChains[ch],algo.MChains[ch].i)
		# oldval = evals(algo.MChains[ch],algo.MChains[ch].i)[1]
		# 1) find all other chains with value +/- x% of chain ch
		close = Int64[]  # vector of indices of "close" chains
		for ch2 in 1:algo["N"]
			if ch != ch2
				e2 = getEval(algo.MChains[ch2],algo.MChains[ch2].i)
				tmp = abs(e2.value - e1.value) / abs(e1.value)
				# tmp = abs(evals(algo.MChains[ch2],algo.MChains[ch2].i)[1] - oldval) / abs(oldval)	# percent deviation
				if tmp < algo.MChains[ch].dist_tol
					push!(close,ch2)
				end
			end
		end
		# 2) with y% probability exchange with a randomly chosen chain from close
		if length(close) > 0
			if rand() < algo.MChains[ch].jump_prob
				ex_with =sample(close)
			#	debug("making an exchange move for chain $ch with chain $ex_with set:$close")
				swapRows!(algo,Pair(ch,ex_with),algo.i)
			end
		end
	end

end

function swapRows!(algo::MAlgoBGP,pair::Pair,i::Int)

	e1 = getEval(algo.MChains[pair.first] ,i)
	e2 = getEval(algo.MChains[pair.second],i)

	# swap
	appendEval!(algo.MChains[pair.first  ],e2)
	appendEval!(algo.MChains[pair.second ],e1)

	# make a note in infos
	algo.MChains[pair.first].infos[i,:exchanged_with] = pair.second
	algo.MChains[pair.second].infos[i,:exchanged_with] = pair.first

end

# get past parameter values to compute new candidates
# we draw new candidates for each chain from a joint normal
# that depends on parameter vectors on ALL chains
# pardf = Dataframe with stacked parameter df for each chain
# function getParamKernel(algo::MAlgoBGP)
function getParamCovariance(algo::MAlgoBGP)

	# select the last algo["past_iterations"] iterations from each chain
	lower_bound_index = maximum([1,algo.MChains[1].i-algo["past_iterations"]])

	# get all params from all chains 
	pardf = parameters(algo.MChains,lower_bound_index:algo.MChains[1].i)
	# |-------|----|------|----------|-----------|
	# | Row # | id | iter | a        | b         |
	# | 1     | 1  | 1    | 0.901685 | 0.983013  |
	# | 2     | 1  | 2    | 0.282522 | 0.763859  |
	# | 3     | 1  | 3    | 0.817773 | 0.268273  |
	# ...

	# compute Var-Cov matrix of parameters_to_sample
	# plus some small random noise
	VV = cov(convert(Array,pardf[:, ps2s_names(algo.m)])) + 0.0001 * Diagonal([1 for i=1:length(ps2s_names(algo.m))])
	return VV

	# setup a MvNormal
	# MVN = MvNormal(VV)
	# return MVN
end


function sample(d::Distributions.MultivariateDistribution,lb::Vector{Float64},ub::Vector{Flaot64})

    # draw as long as all points are in support
    n = 100
    for i in 1:n
        x = rand(d)
        if (x.>=lb) && (x.<=ub)
            return x
        end
    end
    error("no draw in support after $n trials")

end

# function getNewCandidates!(algo::MAlgoBGP,MVN::MvNormal)
function getNewCandidates!(algo::MAlgoBGP,VV::Matrix)

	# update chain by chain
	for ch in 1:algo["N"]

		# TODO
		# getParamKernel could be in here
		# chain.temperature should parameterize MVN somehow (as in their toy example: multiplies the variance)
		# setup a MvNormal


        # Transition Kernel is q(.|theta(t-1)) ~ TruncatedN(theta(t-1), Sigma,lb,ub)
        sig = algo.MChains[ch].shock_sd
        mu  = getLastEval(algo.MChains[ch]).params
        lb = [v[:lb] for (k,v) in algo.m.params_to_sample]
        ub = [v[:ub] for (k,v) in algo.m.params_to_sample]

        newp = sample(MvNormal(mu,sig),lb,ub)


        eval_old = getLastEval(algo.MChains[ch])
        for k in keys(eval_old.params)
            algo.current_param[ch][k] = fitMirror(eval_old.params[k] + shock[k] , 
                                                    algo.m.params_to_sample[k][:lb],
                                                    algo.m.params_to_sample[k][:ub])
        end
		jumpParams!(algo,ch,shock)
		# debug("new params: $(algo.current_param)")
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
function save(algo::MAlgoBGP, filename::AbstractString)
    # step 1, create the file if it does not exist

    ff5 = h5open(filename, "w")

    vals = String[]
    keys = String[]
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

function readAlgoBGP(filename::AbstractString)

    ff5 = h5open(filename, "r")
    keys = HDF5.read(ff5,"algo/opts/keys")
    vals = HDF5.read(ff5,"algo/opts/vals")
    opts = Dict()
    for k in 1:length(keys)
    	opts[keys[k]] = vals[k]
    end

    # each chain has 3 data.frames: parameters, moments and infos
    n = parse(Int,opts["N"])
    params = simpleDataFrameRead(ff5,joinpath("chain","1","parameters"))
    moments = simpleDataFrameRead(ff5,joinpath("chain","1","moments"))
    infos = simpleDataFrameRead(ff5,joinpath("chain","1","infos"))
    if n>1
    	for ich in 2:n
    		params = vcat(params, simpleDataFrameRead(ff5,joinpath("chain","$ich","parameters")))
    		moments = vcat(moments, simpleDataFrameRead(ff5,joinpath("chain","$ich","moments")))
    		infos = vcat(infos, simpleDataFrameRead(ff5,joinpath("chain","$ich","infos")))
    	end
    end
    close(ff5)
    return Dict("opts" => opts, "params"=> params, "moments"=>moments,"infos"=>infos)
end


function show(io::IO,MA::MAlgoBGP)
	print(io,"\n")
	print(io,"BGP Algorithm with $(MA["N"]) chains\n")
	print(io,"============================\n")
	print(io,"\n")
	print(io,"Algorithm\n")
	print(io,"---------\n")
	print(io,"Current iteration: $(MA.i)\n")
	print(io,"Number of params to estimate: $(length(MA.m.params_to_sample))\n")
	print(io,"Number of moments to match: $(length(MA.m.moments))\n")
	print(io,"\n")
	print(io,"Chains\n")
	print(io,"------\n")
	print(io,"Tempering range: [$(MA.MChains[1].tempering),$(MA.MChains[end].tempering)]\n")
	print(io,"Jump probability range: [$(MA.MChains[1].jump_prob),$(MA.MChains[end].jump_prob)]\n")
end
