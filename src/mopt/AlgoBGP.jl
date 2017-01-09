
abstract AbstractChain

# The BGP MCMC Algorithm: Likelihood-Free Parallel Tempering
# ==========================================================
#
# http://link.springer.com/article/10.1007%2Fs11222-012-9328-6
# http://fr.arxiv.org/abs/1108.3423
#
# Baragatti, Grimaud and Pommeret (BGP)
# 
# Approximate Bayesian Computational (ABC) methods (or likelihood-free methods) have appeared in the past fifteen years as useful methods to perform Bayesian analyses when the likelihood is analytically or computationally intractable. Several ABC methods have been proposed: Monte Carlo Markov BGPChains (MCMC) methods have been developped by Marjoramet al. (2003) and by Bortotet al. (2007) for instance, and sequential methods have been proposed among others by Sissonet al. (2007), Beaumont et al. (2009) and Del Moral et al. (2009). Until now, while ABC-MCMC methods remain the reference, sequential ABC methods have appeared to outperforms them (see for example McKinley et al. (2009) or Sisson et al. (2007)). In this paper a new algorithm combining population-based MCMC methods with ABC requirements is proposed, using an analogy with the Parallel Tempering algorithm (Geyer, 1991). Performances are compared with existing ABC algorithms on simulations and on a real example.



###################################
# Start defining BGPChain
###################################

type BGPChain <: AbstractChain
    evals     :: Array{Eval}
    id        :: Int64
    iter      :: Int64
    accepted  :: Array{Bool}
    accept_rate :: Float64
    exchanged :: Array{Int}
    m         :: MProb
    sigma     :: PDiagMat{Float64}
    update_sigma :: Int64   # update sampling vars every x periods
    smpl_iters :: Int64   # max number of trials to get a new parameter from MvNormal that lies within support

    function BGPChain(id::Int,n::Int,m::MProb,sig::Vector{Float64},upd::Int64,smpl_iters::Int)
        @assert length(sig) == length(m.params_to_sample)
        this           = new()
        this.evals     = Array{Eval}(n)
        this.evals[1]  = Eval(m)    # set first eval
        this.accepted  = falses(n)
        this.accept_rate = 0.0
        this.exchanged = zeros(Int,n)
        this.id        = id
        this.iter      = 0
        this.m         = m
        this.sigma     = PDiagMat(sig)
        this.update_sigma = upd
        this.smpl_iters = smpl_iters
        return this
    end
end

function lastAccepted(c::BGPChain)
    if c.iter==1
        return 1
    else
        return find(c.accepted[(c.iter-1):-1:1])[1]
    end
end
getLastEval(c::BGPChain) = c.evals[lastAccepted(c)]
set_sigma!(c::BGPChain,s::Vector{Float64}) = length(s) == length(c.m.params_to_sample) ? c.sigma = PDiagMat(s) : ArgumentError("s has wrong length")
function set_eval!(c::BGPChain,ev::Eval)
    c.evals[c.iter] = deepcopy(ev)
    c.accepted[c.iter] =  ev.accepted 
    return nothing
end
function set_exchanged!(c::BGPChain,i::Int)
    c.exchanged[c.iter] = i
    return nothing
end

function set_acceptRate!(c::BGPChain) 
    c.accept_rate = mean(c.accepted[1:c.iter])
end

function next_eval!(c::BGPChain)
    # generate new parameter vector from last accepted param

    # increment interation
    c.iter += 1

    pp = getNewCandidates(c)

    # evaluate objective 
    ev = evaluateObjective(c.m,pp)

    # accept reject 
    doAcceptReject!(c,ev)

    # save eval on BGPChain 
    set_eval!(c,ev)

end

function doAcceptReject!(c::BGPChain,eval_new::Eval)
    if c.iter == 1
        # accept everything.
        eval_new.prob =1.0
        eval_new.accepted = true
        c.accepted[c.iter] = eval_new.accepted
        set_acceptRate!(c)
    else
        eval_old = getLastEval(c)

        if eval_new.status < 0
            eval_new.prob = 0.0
            eval_new.accepted = false
        else

            eval_new.prob = minimum([1.0,exp( ( eval_old.value - eval_new.value) )])

            if isna(eval_new.prob)
                eval_new.prob = 0.0
                eval_new.accepted = false

            elseif !isfinite(eval_old.value)
                # should never get accepted
                @debug("eval_old is not finite")
                eval_new.prob = 1.0
                eval_new.accepted = true 
            else 
                # status = 1
                if eval_new.prob > rand()
                    eval_new.accepted = true 
                else
                    eval_new.accepted = false
                end
            end 

        end

        c.accepted[c.iter] = eval_new.accepted
        set_acceptRate!(c)

        # update sampling variances every x periods
        # -----------------------------------------

        # update shock variance. want to achieve a long run accpetance rate of 23.4% (See Casella and Berger)
        # and only if you are not BGPChain number 1

        if (c.id>1) && (mod(c.iter,c.update_sigma) == 0)
            too_high = c.accept_rate > 0.234
            if too_high
                @debug("acceptance rate on BGPChain $(c.id) is too high at $(c.accept_rate). increasing variance of each param by 2%.")
                set_sigma!(c,diag(c.sigma) .* 1.02 )
            else
                @debug("acceptance rate on BGPChain $(c.id) is too low at $(c.accept_rate). decreasing variance of each param by 2%.")
                set_sigma!(c,diag(c.sigma) .* 0.98 )
            end
        end
    end
end


"""
    sample(d::Distributions.MultivariateDistribution,lb::Vector{Float64},ub::Vector{Float64},iters::Int)

sample from distribution `d` until all poins are in support
"""
function sample(d::Distributions.MultivariateDistribution,lb::Vector{Float64},ub::Vector{Float64},iters::Int)

    # draw until all points are in support
    n = iters
    for i in 1:n
        x = rand(d)
        if all(x.>=lb) && all(x.<=ub)
            return x
        end
    end
    error("no draw in support after $n trials. increase either opts[smpl_iters] or opts[bound_prob].")
end

function getNewCandidates(c::BGPChain)
    # Transition Kernel is q(.|theta(t-1)) ~ TruncatedN(theta(t-1), Sigma,lb,ub)

    if c.iter==1
        return c.m.initial_value
    else
        ev_old = getLastEval(c)
        mu  = paramd(ev_old) # dict of params
        lb = [v[:lb] for (k,v) in c.m.params_to_sample]
        ub = [v[:ub] for (k,v) in c.m.params_to_sample]

        newp = Dict(zip(collect(keys(mu)),sample(MvNormal(collect(values(mu)),c.sigma),lb,ub,c.smpl_iters)))
        @debug("old param: $(ev_old.params)")
        @debug("new param: $newp")

        return newp
    end

end


###################################
# end BGPChain
###################################



type MAlgoBGP <: MAlgo
    m               :: MProb # an MProb
    opts            :: Dict	# list of options
    i               :: Int 	# iteration
    chains         :: Array{BGPChain} 	# collection of BGPChains: if N==1, length(BGPChains) = 1
  
    function MAlgoBGP(m::MProb,opts=Dict("N"=>3,"maxiter"=>100,"maxtemp"=> 2,"bound_prob"=>0.15,"disttol"=>0.1,"update_sigma"=>10))

        if opts["N"] > 1
    		temps     = linspace(1.0,opts["maxtemp"],opts["N"])
            # initial std dev for each parameter to achieve at least bound_prob on the bounds
            # println("opts=$opts")
            # println("pars = $( m.params_to_sample)")
            init_sd = OrderedDict()
            for (k,v) in m.params_to_sample
                init_sd[k] = MOpt.initvar(v[:lb],(v[:lb]+v[:ub])/2,opts["bound_prob"])
            end
            BGPChains = BGPChain[BGPChain(i,opts["maxiter"],m,collect(values(init_sd)) .* temps[i],opts["update_sigma"],opts["smpl_iters"]) for i in 1:opts["N"]]
          else
            temps     = [1.0]
            init_sd = OrderedDict()
            for (k,v) in m.params_to_sample
                init_sd[k] = MOpt.initvar(v[:lb],(v[:lb]+v[:ub])/2,opts["bound_prob"],opts["smpl_iters"])
            end
            BGPChains = BGPChain[BGPChain(i,opts["maxiter"],m,collect(values(init_sd)),opts["update_sigma"],opts["smpl_iters"])]
        end
	    return new(m,opts,0,BGPChains)
    end
end






# computes new candidate vectors for each BGPChain
# accepts/rejects that vector on each BGPChain, according to some rule
# *) computes N new parameter vectors
# *) applies a criterion to accept/reject any new params
# *) stores the result in BGPChains
function computeNextIteration!( algo::MAlgoBGP )
    # here is the meat of your algorithm:
    # how to go from p(t) to p(t+1) ?

    incrementBGPChainIter!(algo.chains)

    # check algo index is the same on all BGPChains
    for ic in 1:algo["N"]
        @assert algo.i == algo.chains[ic].i
    end

    # TODO 
    # this on each BGPChain
    # START=========================================================
    pmap( x->next_eval!(x), algo.chains ) # this does getNewCandidates, evaluateObjective, doAcceptRecject

    # Part 2) EXCHANGE MOVES only on master
    # ----------------------
    # starting mixing in period 3
    if algo.i>=2 && algo["N"] > 1 
        exchangeMoves!(algo)
    end
end


function exchangeMoves!(algo::MAlgoBGP)
	# for all BGPChains
    dtol = get(algo.opts,"disttol",0.1)
	for ch in 1:algo["N"]
        e1 = getLastEval(algo.chains[ch])
		# 1) find all other BGPChains with value +/- x% of BGPChain ch
		close = Int64[]  # vector of indices of "close" BGPChains
		for ch2 in 1:algo["N"]
			if ch != ch2
				e2 = getLastEval(algo.chains[ch2])
				tmp = abs(e2.value - e1.value) / abs(e1.value)
				# tmp = abs(evals(algo.chains[ch2],algo.chains[ch2].i)[1] - oldval) / abs(oldval)	# percent deviation
				if tmp < dtol 
                    @debug("perc dist $ch and $ch2 is $tmp. label close.")
					push!(close,ch2)
				end
			end
		end
		# 2) with y% probability exchange with a randomly chosen BGPChain from close
		if length(close) > 0
			ex_with = rand(close)
			@debug("making an exchange move for BGPChain $ch with BGPChain $ex_with set:$close")
			swap_ev!(algo,Pair(ch,ex_with),algo.i)
		end
	end

end

function swap_ev!(algo::MAlgoBGP,ch_id::Pair)
    @debug("swapping evs for $(ch_id.first) and $(ch_id.second)")
    c1 = algo.chains[ch_id.first] 
    c2 = algo.chains[ch_id.second]

	e1 = getLastEval(c1)
	e2 = getLastEval(c2)

    # swap 
    set_eval!(c1,e2)
    set_eval!(c2,e1)

	# make a note 
    set_exchanged!(c1,ch_id.second)
    set_exchanged!(c2,ch_id.first)
end






# save algo BGPChains component-wise to HDF5 file
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

	# saving the BGPChains
	for cc in 1:algo["N"]
	    saveBGPChainToHDF5(algo.chains[cc], ff5, "BGPChain/$cc")
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

    # each BGPChain has 3 data.frames: parameters, moments and infos
    n = parse(Int,opts["N"])
    params = simpleDataFrameRead(ff5,joinpath("BGPChain","1","parameters"))
    moments = simpleDataFrameRead(ff5,joinpath("BGPChain","1","moments"))
    infos = simpleDataFrameRead(ff5,joinpath("BGPChain","1","infos"))
    if n>1
    	for ich in 2:n
    		params = vcat(params, simpleDataFrameRead(ff5,joinpath("BGPChain","$ich","parameters")))
    		moments = vcat(moments, simpleDataFrameRead(ff5,joinpath("BGPChain","$ich","moments")))
    		infos = vcat(infos, simpleDataFrameRead(ff5,joinpath("BGPChain","$ich","infos")))
    	end
    end
    close(ff5)
    return Dict("opts" => opts, "params"=> params, "moments"=>moments,"infos"=>infos)
end


function show(io::IO,MA::MAlgoBGP)
	print(io,"\n")
	print(io,"BGP Algorithm with $(MA["N"]) BGPChains\n")
	print(io,"============================\n")
	print(io,"\n")
	print(io,"Algorithm\n")
	print(io,"---------\n")
	print(io,"Current iteration: $(MA.i)\n")
	print(io,"Number of params to estimate: $(length(MA.m.params_to_sample))\n")
	print(io,"Number of moments to match: $(length(MA.m.moments))\n")
	print(io,"\n")
	print(io,"BGPChains\n")
	print(io,"------\n")
	print(io,"Tempering range: [$(MA.chains[1].tempering),$(MA.chains[end].tempering)]\n")
	print(io,"Jump probability range: [$(MA.chains[1].jump_prob),$(MA.chains[end].jump_prob)]\n")
end
