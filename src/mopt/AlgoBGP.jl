
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
type BGPChain <: AbstractChain
	id::Int             # chain id
	i::Int              # current index
	infos        ::DataFrame   # DataFrameionary of arrays(L,1) with eval, ACC and others
	parameters   ::DataFrame   # DataFrameionary of arrays(L,1), 1 for each parameter
	moments      ::DataFrame   # DataFrameionary of DataArrays(L,1), 1 for each moment
	dist_tol     ::Float64 # percentage value distance from current chain i that is considered "close" enough. i.e. close = ((val_i - val_j)/val_j < dist_tol )
	jump_prob    ::Float64 # probability of swapping with "close" chain

	params_nms   ::Array{Symbol,1}	# names of parameters (i.e. exclusive of "id" or "iter", etc)
	moments_nms  ::Array{Symbol,1}	# names of moments
	params2s_nms ::Array{Symbol,1}  # DataFrame names of parameters to sample 

	tempering  ::Float64 # tempering in update probability
	shock_sd   ::Float64 # sd of shock to 

	function BGPChain(id,MProb,L,temp,shock,dist_tol,jump_prob)
		infos      = DataFrame(chain_id = [id for i=1:L], iter=1:L, evals = DataArray(zeros(Float64,L)), accept = zeros(Bool,L), status = zeros(Int,L), exchanged_with=zeros(Int,L),prob=zeros(Float64,L),perc_new_old=zeros(Float64,L),accept_rate=zeros(Float64,L),shock_sd = [shock;zeros(Float64,L-1)],eval_time=zeros(Float64,L),tempering=zeros(Float64,L))
		parameters = DataFrame(chain_id = [id for i=1:L], iter=1:L)
        moments    = DataFrame(chain_id = [id for i=1:L], iter=1:L)
		par_nms    = sort(Symbol[ symbol(x) for x in ps_names(MProb) ])
		par2s_nms  = Symbol[ symbol(x) for x in ps2s_names(MProb) ]
		mom_nms    = sort(Symbol[ symbol(x) for x in ms_names(MProb) ])
        for i in par2s_nms
            parameters[i] = DataArray(zeros(L))
        end
        for i in mom_nms
            moments[i] = DataArray(zeros(L))
        end
		return new(id,0,infos,parameters,moments,dist_tol,jump_prob,par_nms,mom_nms,par2s_nms,temp,shock)
    end
end

type BGPChains
	MChains :: Array{BGPChain,1}
end

type MAlgoBGP <: MAlgo
    m               :: MProb # an MProb
    opts            :: Dict	# list of options
    i               :: Int 	# iteration
    current_param   :: Array{Dict,1}  # current param value: one Dict for each chain
    MChains         :: Array{BGPChain,1} 	# collection of Chains: if N==1, length(chains) = 1
  
    function MAlgoBGP(m::MProb,opts=Dict("N"=>3,"min_shock_sd"=>0.1,"max_shock_sd"=>1.0,"maxiter"=>100,"maxtemp"=> 100))

        if opts["N"] > 1
    		temps     = linspace(1.0,opts["maxtemp"],opts["N"])
    		shocksd   = linspace(opts["min_shock_sd"],opts["max_shock_sd"],opts["N"])
    		disttol   = linspace(opts["min_disttol"],opts["max_disttol"],opts["N"])
    		jump_prob = linspace(opts["min_jump_prob"],opts["max_jump_prob"],opts["N"])
    	  	chains    = [BGPChain(i,m,opts["maxiter"],temps[i],shocksd[i],disttol[i],jump_prob[i]) for i=1:opts["N"] ]
          else
            temps     = [1.0]
            shocksd   = [opts["min_shock_sd"]]
            disttol   = [opts["min_disttol"]]
            jump_prob = [opts["min_jump_prob"]]
            chains    = [BGPChain(i,m,opts["maxiter"],temps[i],shocksd[i],disttol[i],jump_prob[i]) for i=1:opts["N"] ]
        end
	  	# current param values
	  	cpar = [ deepcopy(m.initial_value) for i=1:opts["N"] ] 

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

	# Part 2) EXCHANGE MOVES 
	# ----------------------
	# starting mixing in period 3
	if algo.i>=2 && algo["N"] > 1 
		exchangeMoves!(algo)
	end
end

# notice: higher tempering draws candiates further spread out,
# but accepts lower function values with lower probability
function doAcceptRecject!(algo::MAlgoBGP,v::Array)
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
            # -------------------------

            # update average acceptance rate: moving average
		    chain.infos[algo.i,:accept_rate]   = 0.9 * chain.infos[algo.i-1,:accept_rate] + 0.1 * ACC

            # update shock variance. want to achieve a long run accpetance rate of 23.4% (See Casella and Berger)
            accept_too_high = chain.infos[algo.i,:accept_rate]>0.234
            if accept_too_high
                chain.shock_sd *= 1.05  # increase variance by 5% => will accept less
            else # too low
                chain.shock_sd *= 0.95  # decrease variance by 5% => will accept more
            end
		    chain.infos[algo.i,:shock_sd]      = chain.shock_sd
		    chain.infos[algo.i,:perc_new_old] = (eval_new.value - eval_old.value) / abs(eval_old.value)
		    # debug("ACCEPTED: $ACC")
		    # debug("old value: $(eval_old.value)")
		    # debug("new value: $(eval_new.value)")
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
		shock_list = [ (algo.m.params_to_sample[p][:ub] - algo.m.params_to_sample[p][:lb]) for p in ps2s_names(algo) ] ./ (1.96 * 2 *diag(VV2))
		shock_ub  = minimum(   shock_list  )
		algo.MChains[ch].shock_sd  = min(algo.MChains[ch].shock_sd , shock_ub)

		# shock parameters on chain index ch
		shock = rand(MVN) * algo.MChains[ch].shock_sd
		shock = Dict(zip(ps2s_names(algo) , shock))

		# debug("shock to parameters on chain $ch :")
		# debug("shock = $shock")

		# debug("current params: $(getLastEval(algo.MChains[ch]).params)")

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
