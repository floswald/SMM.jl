
abstract type AbstractChain end





###################################
# Start defining BGPChain
###################################




"""
# `BGPChain`

MCMC Chain storage for BGP algorithm. This is the main datatype for the implementation of Baragatti, Grimaud and Pommeret (BGP) in [Likelihood-free parallel tempring](http://arxiv.org/abs/1108.3423)

## Fields

* `evals`: Array of `Eval`s
* `best_id`: index of best `eval.value` so far
* `best_val`: best eval.value so far
* `curr_val` : current value
* `probs_acc`: vector of probabilities with which to accept current value
* `id`: Chain identifier
* `iter`: current iteration
* `accepted`: `Array{Bool}` of `length(evals)`
* `accept_rate`: current acceptance rate
* `acc_tuner`: Acceptance tuner. `acc_tuner > 1` means to be more restrictive: params that yield a *worse* function value are *less likely* to get accepted, the higher `acc_tuner` is.
* `exchanged`: `Array{Int}` of `length(evals)` with index of chain that was exchanged with
* `m`: `MProb`
* `sigma`: `Float64` shock variance
* `sigma_update_steps`:  update sampling vars every `sigma_update_steps` iterations. setting `sigma_update_steps > maxiter` means to never update the variances.
* `sigma_adjust_by`: adjust sampling vars by `sigma_adjust_by` percent up or down
* `smpl_iters`: max number of trials to get a new parameter from MvNormal that lies within support
* `min_improve`: minimally required improvement in chain `j` over chain `i` for an exchange move `j->i` to talk place.
* `batches`: in the proposal function update the parameter vector in batches. [default: update entire param vector]

"""
mutable struct BGPChain <: AbstractChain
    evals     :: Array{Eval}
    best_id   :: Vector{Int}   # index of best eval.value so far
    best_val  :: Vector{Float64}   # best eval.value so far
    curr_val  :: Vector{Float64}   # current value
    probs_acc :: Vector{Float64}    # vector of probabilities with which to accept
    id        :: Int64
    iter      :: Int64
    accepted  :: Array{Bool}
    accept_rate :: Float64
    acc_tuner :: Float64
    exchanged :: Array{Int}
    m         :: MProb
    sigma     :: Float64
    sigma_update_steps :: Int64   # update sampling vars every sigma_update_steps iterations
    sigma_adjust_by :: Float64   # adjust sampling vars by sigma_adjust_by percent up or down
    smpl_iters :: Int64   # max number of trials to get a new parameter from MvNormal that lies within support
    min_improve  :: Float64  
    batches  :: Vector{UnitRange{Int}}  # vector of indices to update together.

    """
        BGPChain(id::Int=1,n::Int=10;
            m::MProb=MProb(),sig::Float64=0.5,upd::Int64=10,upd_by::Float64=0.01,smpl_iters::Int=1000,
            min_improve::Float64=10.0,acc_tuner::Float64=2.0,batch_size=1)

    Constructor of a BGPChain. Keyword args:
        * `acc_tuner`: Acceptance tuner. `acc_tuner > 1` means to be more restrictive: params that yield a *worse* function value are *less likely* to get accepted, the higher `acc_tuner` is.
        * `exchanged`: `Array{Int}` of `length(evals)` with index of chain that was exchanged with
        * `m`: `MProb`
        * `sig`: `Float64` shock variance
        * `upd`:  update sampling vars every `upd` iterations
        * `upd_by`: adjust sampling vars by `upd_by` percent up or down
        * `smpl_iters`: max number of trials to get a new parameter from MvNormal that lies within support
        * `min_improve`: minimally required improvement in chain `j` over chain `i` for an exchange move `j->i` to talk place.
        * `batch_size`: size of batches in which to update parameter vector.
    """
    function BGPChain(id::Int=1,n::Int=10;m::MProb=MProb(),sig::Float64=0.5,upd::Int=10,upd_by::Float64=0.01,smpl_iters::Int=1000,min_improve::Float64=10.0,acc_tuner::Float64=2.0,batch_size=1)
        np = length(m.params_to_sample)
        this           = new()
        this.evals     = Array{Eval}(undef,n)
        this.best_val  = ones(n) * Inf
        this.best_id   = -ones(Int,n)
        this.curr_val  = ones(n) * Inf
        this.probs_acc = rand(n)
        this.evals[1]  = Eval(m)    # set first eval
        this.accepted  = falses(n)
        this.accept_rate = 0.0
        this.acc_tuner = acc_tuner
        this.exchanged = zeros(Int,n)
        this.id        = id
        this.iter      = 0
        this.m         = m
        # how many bundles + rest
        nb, rest = divrem(np,batch_size)
        this.sigma = sig 
        this.batches = UnitRange{Int}[]
        i = 1
        for ib in 1:nb
            j = (ib==nb && rest > 0) ? length(sig) :  i + batch_size - 1
            push!(this.batches,i:j)
            i = j + 1
        end
        this.sigma_update_steps = upd
        this.sigma_adjust_by = upd_by
        this.smpl_iters = smpl_iters
        this.min_improve = min_improve
        return this
    end
end

"""
    allAccepted(c::BGPChain)

Get all accepted `Eval`s from a chain
"""
allAccepted(c::BGPChain) = c.evals[c.accepted]

# return a dict of param values as arrays
function params(c::BGPChain;accepted_only=true)
    if accepted_only
        e = allAccepted(c)
    else
        e = c.evals
    end
    d = Dict{Symbol,Vector{Float64}}()
    for k in keys(e[1].params)
        d[k] = Float64[e[i].params[k] for i in 1:length(e)]
    end
    return d
end

"""
    history(c::BGPChain)

Returns a `DataFrame` with a history of the chain.
"""
function history(c::BGPChain)
    N = length(c.evals)
    cols = Any[]
    # d = DataFrame([Int64,Float64,Bool,Int64],[:iter,:value,:accepted,:prob],N)
    d = DataFrame()
    d.iter = collect(1:c.iter)
    d[!,:exchanged] = c.exchanged
    d[!,:accepted] = c.accepted
    d[!,:best_val] = c.best_val
    d[!,:curr_val] = c.curr_val
    d[!,:best_id] = c.best_id
    # get fields from evals
    nms = [:value,:prob]
    for n in nms
        d[!,n] = eltype(getfield(c.evals[1],n))[getfield(c.evals[i],n) for i in 1:N]
    end
    # get fields from evals.params
    for (k,v) in c.evals[1].params
        d[!,k] = eltype(v)[c.evals[i].params[k] for i in 1:N]
    end

    return d[!,[:iter,:value,:accepted,:curr_val, :best_val, :prob, :exchanged,collect(keys(c.evals[1].params))...]]
end

"""
    best(c::BGPChain) -> (val,idx)

Returns the smallest value and index stored of the chain.
"""
best(c::BGPChain) = findmin([c.evals[i].value for i in 1:length(c.evals)])

"""
    mean(c::BGPChain)

Returns the mean of all parameter values stored on the chain.
"""
mean(c::BGPChain) = Dict(k => mean(v) for (k,v) in params(c))

"""
    median(c::BGPChain)

Returns the median of all parameter values stored on the chain.
"""
median(c::BGPChain) = Dict(k => median(v) for (k,v) in params(c))

"""
    CI(c::BGPChain;level=0.95)

Confidence interval on parameters
"""
CI(c::BGPChain;level=0.95) = Dict(k => quantile(v,[(1-level)/2, 1-(1-level)/2]) for (k,v) in params(c))



"""
    summary(c::BGPChain)

Returns a summary of the chain. Condensed [`history`](@ref)
"""
function summary(c::BGPChain)
    ex_with = c.exchanged[c.exchanged .!= 0]
    if length(ex_with) == 1
        ex_with  = [ex_with]
    end
    d = DataFrame(id =c.id, acc_rate = c.accept_rate,perc_exchanged=100*sum(c.exchanged .!= 0)/length(c.exchanged),
        exchanged_most_with=length(ex_with)>0 ? mode(ex_with) : 0,
        best_val=c.best_val[end])
    return d
end


function lastAccepted(c::BGPChain)
    if c.iter==1
        return 1
    else
        return findlast(c.accepted[1:(c.iter)])
    end
end
getIterEval(c::BGPChain,i::Int) = c.evals[i]
getLastAccepted(c::BGPChain) = c.evals[lastAccepted(c)]
# set_sigma!(c::BGPChain,s::Vector{Float64}) = length(s) == length(c.m.params_to_sample) ? c.sigma = PDiagMat(s) : ArgumentError("s has wrong length")
set_sigma!(c::BGPChain,s::Float64) = c.sigma = s
function set_eval!(c::BGPChain,ev::Eval)
    c.evals[c.iter] = deepcopy(ev)
    c.accepted[c.iter] =  ev.accepted
    # set best value
    if c.iter == 1
        c.best_val[c.iter] = ev.value
        c.curr_val[c.iter] = ev.value
        c.best_id[c.iter] = c.iter

    else
        if ev.accepted
            c.curr_val[c.iter] = ev.value
        else
            c.curr_val[c.iter] = c.curr_val[c.iter-1]
        end
        if (ev.value < c.best_val[c.iter-1])
            c.best_val[c.iter] = ev.value
            c.best_id[c.iter]  = c.iter
        else
            # otherwise, keep best and current from last iteration
            c.best_val[c.iter] = c.best_val[c.iter-1]
            c.best_id[c.iter]  = c.best_id[c.iter-1]
        end
    end
    return nothing
end
function set_exchanged!(c::BGPChain,i::Int)
    c.exchanged[c.iter] = i
    return nothing
end


"set acceptance rate on chain. considers only iterations where no exchange happened."
function set_acceptRate!(c::BGPChain)
    noex = c.exchanged[1:c.iter] .== 0
    acc = c.accepted[1:c.iter]
    c.accept_rate = mean(acc[noex])
end


"""
    next_eval(c::BGPChain)

Computes the next `Eval` for chain `c`:

1. Get last accepted param 
2. get a new param via [`proposal`](@ref)
3. [`evaluateObjective`](@ref)
4. Accept or Reject the new value via [`doAcceptReject!`](@ref)
5. Store `Eval` on chain `c`.

"""
function next_eval(c::BGPChain)
    # generate new parameter vector from last accepted param

    # increment interation
    c.iter += 1
    @debug "iteration = $(c.iter)"

    # returns an OrderedDict
    pp = proposal(c)

    # evaluate objective
    ev = evaluateObjective(c.m,pp)


    # accept reject
    doAcceptReject!(c,ev)

    # save eval on BGPChain
    set_eval!(c,ev)

    return c

end


"""
    doAcceptReject!(c::BGPChain,eval_new::Eval)

Perform a Metropolis-Hastings accept-reject operation on the latest `Eval` and update the sampling variance, if so desired (set via `sigma_update_steps` in [`BGPChain`](@ref) constructor.)
"""
function doAcceptReject!(c::BGPChain,eval_new::Eval)
            @debug "doAcceptReject!"
    if c.iter == 1
        # accept everything.
        eval_new.prob =1.0
        eval_new.accepted = true
        eval_new.status = 1
        c.accepted[c.iter] = eval_new.accepted
        set_acceptRate!(c)
    else
        eval_old = getLastAccepted(c)

        if eval_new.status < 0
            eval_new.prob = 0.0
            eval_new.accepted = false
        else

            eval_new.value >= 0 || error("AlgoBGP assumes that your objective function returns a non-negative number.")
            # this forumulation: old - new
            # because we are MINIMIZING the value of the objective function
            eval_new.prob = minimum([1.0,exp( c.acc_tuner * ( eval_old.value - eval_new.value) )]) #* (eval_new.value < )
            @debug "eval_new.value = $(eval_new.value)"
            @debug "eval_old.value = $(eval_old.value)"
            @debug "eval_new.prob = $(round(eval_new.prob,digits = 2))"
            @debug "c.probs_acc[c.iter] = $(round(c.probs_acc[c.iter],digits = 2))"

            if !isfinite(eval_new.prob)
                eval_new.prob = 0.0
                eval_new.accepted = false
                eval_new.status = -1

            elseif !isfinite(eval_old.value)
                # should never have gotten accepted
                @debug "eval_old is not finite"
                eval_new.prob = 1.0
                eval_new.accepted = true
            else
                # status = 1
                eval_new.status = 1
                if eval_new.prob > c.probs_acc[c.iter]
                    eval_new.accepted = true
                else
                    eval_new.accepted = false
                end
            end
            @debug "eval_new.accepted = $(eval_new.accepted)"

        end

        c.accepted[c.iter] = eval_new.accepted
        set_acceptRate!(c)

        # update sampling variances every x periods
        # -----------------------------------------

        # update shock variance. want to achieve a long run accpetance rate of 23.4% (See Casella and Berger)

        if mod(c.iter,c.sigma_update_steps) == 0
            too_high = c.accept_rate > 0.234
            if too_high
                @debug "acceptance rate on BGPChain $(c.id) is too high at $(c.accept_rate). increasing variance of each param by $(100* c.sigma_adjust_by)%."
                set_sigma!(c,c.sigma .* (1.0+c.sigma_adjust_by) )
            else
                @debug "acceptance rate on BGPChain $(c.id) is too low at $(c.accept_rate). decreasing variance of each param by $(100* c.sigma_adjust_by)%."
                set_sigma!(c,c.sigma .* (1.0-c.sigma_adjust_by) )
            end
        end
    end
end


"""
Truncated multivariate normal distribution.

From https://github.com/JuliaStats/Distributions.jl/issues/480#issuecomment-1316525345
"""
function truncated_mv_normal(mu::Vector{Float64}, sigma::Float64, lb::Float64, ub::Float64)
    Product([truncated(Normal(m,sigma), lb, ub) for m in mu])
end

"""
    proposal(c::BGPChain)

Gaussian Transition Kernel centered on current parameter value. 

1. Map all ``k`` parameters into ``\\mu \\in [0,1]^k``.
2. update all parameters by sampling from `MvNormal`, ``N(\\mu,\\sigma)``, where ``sigma`` is `c.sigma` until all params are in ``[0,1]^k``
3. Map ``[0,1]^k`` back to original parameter spaces.

"""
function proposal(c::BGPChain)

    if c.iter==1
        return c.m.initial_value
    else
        ev_old = getLastAccepted(c)
        mu  = paramd(ev_old) # dict of params
        lb = [v[:lb] for (k,v) in c.m.params_to_sample]
        ub = [v[:ub] for (k,v) in c.m.params_to_sample]

        # map into [0,1]
        # (x-a)/(b-a) = z \in [0,1]
        mu01 = mapto_01(mu,lb,ub) 
        
        # Transition Kernel is q(.|theta(t-1)) ~ TruncatedN(theta(t-1), Sigma,lb,ub)

        # if there is only one batch of params
        if length(c.batches) == 1
            d = truncated_mv_normal(mu01, c.sigma, 0.0, 1.0)
            pp = mysample(d,0.0,1.0,c.smpl_iters)
        else
            # do it in batches of params
            pp = zero(mu01)
            for (sig_ix,i) in enumerate(c.batches)
                try
                    d = truncated_mv_normal([mu01[i]], c.sigma, 0.0, 1.0)
                    pp[i] = mysample(d,0.0,1.0,c.smpl_iters)
                catch err
                    @error "caught exception $err. this is param index $sig_ix, mean = $(mu01[i]), sigma $(c.sigma), lb,ub = $((0,1))"
                end
            end
        end
        # map [0,1] -> [a,b]
        # z*(b-a) + a = x \in [a,b]
        newp = OrderedDict(zip(collect(keys(mu)),mapto_ab(pp,lb,ub)))

        @debug "iteration $(c.iter)"
        @debug "old param: $(ev_old.params)"
        @debug "new param: $newp"
        for (k,v) in newp
            @debug "step for $k = $(v-ev_old.params[k])"
        end

        # flat kernel: random choice in each dimension.
        # newp = Dict(zip(collect(keys(mu)),rand(length(lb)) .* (ub .- lb)))

        return newp
    end

end


###################################
# end BGPChain
###################################


"""
# MAlgoBGP: BGP MCMC Algorithm

This implements the [BGP MCMC Algorithm Likelihood-Free Parallel Tempering](http://fr.arxiv.org/abs/1108.3423) by Baragatti, Grimaud and Pommeret (BGP):

> Approximate Bayesian Computational (ABC) methods (or likelihood-free methods) have appeared in the past fifteen years as useful methods to perform Bayesian analyses when the likelihood is analytically or computationally intractable. Several ABC methods have been proposed: Monte Carlo Markov BGPChains (MCMC) methods have been developped by Marjoramet al. (2003) and by Bortotet al. (2007) for instance, and sequential methods have been proposed among others by Sissonet al. (2007), Beaumont et al. (2009) and Del Moral et al. (2009). Until now, while ABC-MCMC methods remain the reference, sequential ABC methods have appeared to outperforms them (see for example McKinley et al. (2009) or Sisson et al. (2007)). In this paper a new algorithm combining population-based MCMC methods with ABC requirements is proposed, using an analogy with the Parallel Tempering algorithm (Geyer, 1991). Performances are compared with existing ABC algorithms on simulations and on a real example.


## Fields

* `m`: [`MProb`](@ref)
* `opts`: a `Dict` of options
* `i`: current iteration
* `chains`: An array of [`BGPChain`](@ref)
* `anim`: `Plots.Animation`
* `dist_fun`: function to measure distance between one evaluation and the next.

"""
mutable struct MAlgoBGP <: MAlgo
    m               :: MProb # an MProb
    opts            :: Dict	# list of options
    i               :: Int 	# iteration
    chains         :: Array{BGPChain} 	# collection of BGPChains: if N==1, length(BGPChains) = 1
    anim           :: Plots.Animation
    dist_fun   :: Function

    function MAlgoBGP(m::MProb,opts=Dict("N"=>3,"maxiter"=>100,"maxtemp"=> 2,"sigma"=>0.05,"sigma_update_steps"=>10,"sigma_adjust_by"=>0.01,"smpl_iters"=>1000,"parallel"=>false,"min_improve"=>[0.0 for i in 1:3],"acc_tuners"=>[2.0 for i in 1:3]))

        if opts["N"] > 1
    		temps     = range(1.0,stop=opts["maxtemp"],length=opts["N"])
            # initial std dev for each parameter to achieve at least bound_prob on the bounds
            # println("opts=$opts")
            # println("pars = $( m.params_to_sample)")

            # choose inital sd for each parameter p
            # such that Pr( x \in [init-b,init+b]) = 0.975
            # where b = (p[:ub]-p[:lb])*opts["coverage"] i.e. the fraction of the search interval you want to search around the initial value
            BGPChains = BGPChain[BGPChain(i,opts["maxiter"],
                m = m,
                sig = get(opts,"sigma",0.05) .* temps[i],
                upd = get(opts,"sigma_update_steps",10),
                upd_by = get(opts,"sigma_adjust_by",0.01),
                smpl_iters = get(opts,"smpl_iters",1000),
                min_improve = get(opts,"min_improve",[0.5 for j in 1:opts["N"]])[i],
                acc_tuner = get(opts,"acc_tuners",[2.0 for j in 1:opts["N"]])[i],
                batch_size = get(opts,"batch_size",length(m.params_to_sample))) for i in 1:opts["N"]]
        else
            # println(init_sd)
            BGPChains = BGPChain[BGPChain(1,opts["maxiter"],
                m = m,
                sig = get(opts,"sigma",0.05),
                upd = get(opts,"sigma_update_steps",10),
                upd_by = get(opts,"sigma_adjust_by",0.01),
                smpl_iters = get(opts,"smpl_iters",1000),
                min_improve = get(opts,"min_improve",[0.5 for j in 1:opts["N"]])[i],
                acc_tuner = get(opts,"acc_tuners",[2.0 for j in 1:opts["N"]])[i],
                batch_size = get(opts,"batch_size",length(m.params_to_sample))) for i in 1:opts["N"]]
        end
	    return new(m,opts,0,BGPChains, Animation(),get(opts,"dist_fun",-))
    end
end

function summary(m::MAlgoBGP)
    s = map(summary,m.chains)
    df = s[1]
    if length(s) > 1
        for i in 2:length(s)
            df = vcat(df,s[i])
        end
    end
    return df
end





# return current param spaces on algo
cur_param(m::MAlgoBGP) = iter_param(m,m.i)


# return param spaces on algo at iter
function iter_param(m::MAlgoBGP,iter::Int)
    r = Dict()
    for ic in 1:length(m.chains)
        if m.i == 0
            r[ic] = Dict(:mu => m.m.initial_value,:sigma => m.chains[ic].sigma)
        else
            ev_old = getIterEval(m.chains[ic],iter)
            r[ic] = Dict(:mu => paramd(ev_old),:sigma => m.chains[ic].sigma)
        end
    end
    r

end



"""
    computeNextIteration!( algo::MAlgoBGP )

computes new candidate vectors for each [`BGPChain`](@ref)
accepts/rejects that vector on each BGPChain, according to some rule. The evaluation objective functions is performed in parallel, is so desired.

1. On each chain `c`:
    * computes new parameter vectors
    * applies a criterion to accept/reject any new params
    * stores the result in BGPChains
2. Calls [`exchangeMoves!`](@ref) to swap chains
"""
function computeNextIteration!( algo::MAlgoBGP )
    # here is the meat of your algorithm:
    # how to go from p(t) to p(t+1) ?

    # incrementBGPChainIter!(algo.chains)


    # TODO
    # this is probably inefficeint
    # ideally, would only pmap evaluateObjective, to avoid large data transfers to each worker (now we're transferring each chain back and forth to each worker.)

    # if get(algo.opts, "parallel", false)
        cs = pmap( x->next_eval(x), algo.chains ) # this does proposal, evaluateObjective, doAcceptRecject
    # else
        # for i in algo.chains
        #     @debug(logger," ")
        #     @debug(logger," ")
        #     @debug(logger,"debugging chain id $(i.id)")
        #     next_eval!(i)
        # end
        # cs = map( x->next_eval(x), algo.chains ) # this does proposal, evaluateObjective, doAcceptRecject
    # end
    # reorder and insert into algo
    for i in 1:algo.opts["N"]
        algo.chains[i] = cs[map(x->getfield(x,:id) == i,cs)][1]
        @assert algo.chains[i].id == i
    end
    if get(algo.opts, "animate", false)
        p = plot(algo,1);
        frame(algo.anim)
    end
    # p = plot(algo,1)
    # display(p)
    # sleep(.1)

    # check algo index is the same on all BGPChains
    for ic in 1:algo["N"]
        @assert algo.i == algo.chains[ic].iter
    end

    # Part 2) EXCHANGE MOVES only on master
    # ----------------------
    # starting mixing in period 3
    if algo.i>=2 && algo["N"] > 1
        exchangeMoves!(algo)
    end
end

"""
    exchangeMoves!(algo::MAlgoBGP)

Exchange chain `i` and `j` with if `dist_fun(evi.value,evj.value)` is greate than a threshold value `c.min_improve`. Commonly, this means that we only exchange if `j` is better by *at least* `c.min_improve`.
"""
function exchangeMoves!(algo::MAlgoBGP)

    # i is new index
    # j is old index

    # algo["N"] exchange moves are proposed
    props = [(i,j) for i in 1:algo["N"], j in 1:algo["N"] if (i<j)]
    # N pairs of chains are chosen uniformly in all possibel pairs with replacement
    samples = algo["N"] < 3 ? algo["N"]-1 : algo["N"]
    pairs = sample(props,samples,replace=false)

    # @debug(logger,"")
    # @debug(logger,"exchangeMoves: proposing pairs")
    # @debug(logger,"$pairs")

    for p in pairs
        i,j = p
        evi = getLastAccepted(algo.chains[i])
        evj = getLastAccepted(algo.chains[j])
        # my version
        # if rand() < algo["mixprob"]
            # if (evj.value < evi.value)  # if j's value is better than i's
            #     @debug(logger,"$j better than $i")
            #     # @debug(logger,"$(abs(j.value)) < $(algo.chains[p[1]].min_improve)")
            #     # swap_ev!(algo,p)
            #     set_ev_i2j!(algo,i,j)
            # else
            #     @debug(logger,"$i better than $j")
            #     set_ev_i2j!(algo,j,i)
            # end
        # end

        # BGP version
        # exchange i with j if rho(S(z_j),S(data)) < epsilon_i
        # @debug(logger,"Exchanging $i with $j? Distance is $(algo.dist_fun(evj.value, evi.value))")
        # @debug(logger,"Exchange: $(algo.dist_fun(evj.value, evi.value)  < algo.chains[i].min_improve)")
        # println("Exchanging $i with $j? Distance is $(algo.dist_fun(evj.value, evi.value))")
        # println("Exchange: $(algo.dist_fun(evj.value, evi.value)  > algo["min_improve"][i])")
        # this formulation assumes that evi.value > 0 always, for all i.
        # swap for sure if there is an improvement, i.e. algo.dist_fun(evj.value, evi.value) > 0
        # swap even if there is a deterioration, but only up to threshold min_improve[i]
        if algo.dist_fun(evi.value, evj.value) > algo.chains[i].min_improve
            swap_ev_ij!(algo,i,j)
        end
    end

	# for ch in 1:algo["N"]
 #        e1 = getLastAccepted(algo.chains[ch])
	# 	# 1) find all other BGPChains with value +/- x% of BGPChain ch
	# 	close = Int64[]  # vector of indices of "close" BGPChains
	# 	for ch2 in 1:algo["N"]
	# 		if ch != ch2
	# 			e2 = getLastAccepted(algo.chains[ch2])
	# 			tmp = abs(e2.value - e1.value) / abs(e1.value)
	# 			# tmp = abs(evals(algo.chains[ch2],algo.chains[ch2].i)[1] - oldval) / abs(oldval)	# percent deviation
	# 			if tmp < dtol
 #                    @debug(logger,"perc dist $ch and $ch2 is $tmp. will label that `close`.")
	# 				push!(close,ch2)
	# 			end
	# 		end
	# 	end
	# 	# 2) with y% probability exchange with a randomly chosen BGPChain from close
	# 	if length(close) > 0
	# 		ex_with = rand(close)
	# 		@debug(logger,"making an exchange move for BGPChain $ch with BGPChain $ex_with set: $close")
	# 		swap_ev!(algo,Pair(ch,ex_with))
	# 	end
	# end

end

function set_ev_i2j!(algo::MAlgoBGP,i::Int,j::Int)
    @debug "setting ev of $i to ev of $j"
    ci = algo.chains[i]
    cj = algo.chains[j]

    ei = getLastAccepted(ci)
    ej = getLastAccepted(cj)

    # set ei -> ej
    set_eval!(ci,ej)

    # make a note
    set_exchanged!(ci,j)
end

"replace the current `Eval` of chain ``i`` with the one of chain ``j``"
function swap_ev_ij!(algo::MAlgoBGP,i::Int,j::Int)
    @debug "swapping ev of $i with ev of $j"
    ci = algo.chains[i]
    cj = algo.chains[j]

    ei = getLastAccepted(ci)
    ej = getLastAccepted(cj)

    # set ei -> ej
    set_eval!(ci,ej)
    set_eval!(cj,ei)

    # make a note
    set_exchanged!(ci,j)
    set_exchanged!(cj,i)
end



"""
  extendBGPChain!(chain::BGPChain, algo::MAlgoBGP, extraIter::Int64)

Starting from an existing [`MAlgoBGP`](@ref), allow for additional iterations
by extending a specific chain. This function is used to restart a previous estimation run via [`restart!`](@ref)
"""
function extendBGPChain!(chain::BGPChain, algo::MAlgoBGP, extraIter::Int64)

  initialIter = algo.i
  finalIter = algo.i + extraIter

  # stores the original chain:
  #---------------------------
  copyOriginalChain = deepcopy(chain)


  # I have to change the following fields:
  #---------------------------------------
  # 1. Change the size:
  #--------------------
  chain.evals    = Array{Eval}(undef,finalIter)
  chain.best_val  = ones(finalIter) * Inf
  chain.best_id   = -ones(Int, finalIter)
  chain.curr_val  = ones(finalIter) * Inf
  chain.probs_acc = rand(finalIter)
  chain.accepted  = falses(finalIter)
  chain.exchanged = zeros(Int,finalIter)


  # 2. Push the previous values in:
  #--------------------------------
  for iterNumber = 1:initialIter
    chain.evals[iterNumber] = copyOriginalChain.evals[iterNumber]
    chain.best_val[iterNumber] = copyOriginalChain.best_val[iterNumber]
    chain.best_id[iterNumber] =  copyOriginalChain.best_id[iterNumber]
    chain.curr_val[iterNumber] = copyOriginalChain.curr_val[iterNumber]
    chain.probs_acc[iterNumber] =  copyOriginalChain.probs_acc[iterNumber]
    chain.accepted[iterNumber] = copyOriginalChain.accepted[iterNumber]
    chain.exchanged[iterNumber] = copyOriginalChain.exchanged[iterNumber]
  end



end

"""
  restart!(algo::MAlgoBGP, extraIter::Int64)

Starting from an existing AlgoBGP, restart the optimization from where it
stopped. Add `extraIter` additional steps to the optimization process.
"""
function restart!(algo::MAlgoBGP, extraIter::Int64)

  @info "Restarting estimation loop with $(extraIter) iterations."
  @info "Current best value on chain 1 before restarting $(SMM.summary(algo)[:best_val][1])"
  t0 = time()

  # Minus 1, to follow the SMM convention
  initialIter = algo.i
  finalIter = initialIter + extraIter

  # Extend algo.chain[].evals
  #---------------------------
  # Loop over chains
  for chainNumber = 1:algo.opts["N"]
    extendBGPChain!(algo.chains[chainNumber], algo, extraIter)
  end

  # To follow the conventions in SMM:
  #----------------------------------------
  for ic in 1:algo["N"]
      algo.chains[ic].iter =   algo.i - 1
   end

  #change maxiter in the dictionary storing options
  #------------------------------------------------
  @debug "Setting algo.opts[\"maxiter\"] = $(finalIter)"
  algo.opts["maxiter"] = finalIter

  # do iterations, starting at initialIter
  # and not at i=1, as in run!
    @info "restarting estimation"
    @showprogress for i in initialIter:finalIter

    algo.i = i


    # try
      computeNextIteration!( algo )

      # If the user stops the execution (control + C), save algo in a file
      # with a special name.
      # I leave commented for the moment
      #-------------------------------------------------------------------
      # catch x
      #   @warn(logger, "Error = ", x)
      #
      #   if isa(x, InterruptException)
      #     @warn(logger, "User interupted the execution")
      #     @warn(logger, "Saving the algorithm to disk.")
      #     save(algo, "InterruptedAlgorithm")
      #   end
      # end

      # Save the process every $(save_frequency) iterations:
      #----------------------------------------------------
      # save at certain frequency
    	if haskey(algo.opts,"save_frequency") == true
        # if the user provided a filename in the options dictionary
    		if haskey(algo.opts,"filename") == true
				if mod(i,algo.opts["save_frequency"]) == 0
					save(algo,algo.opts["filename"])
					@info(logger,"saved data at iteration $i")
				end
            end
    	end
        t1 = round((time()-t0)/60,digits = 1)
    	algo.opts["time"] = t1
    	if haskey(algo.opts,"filename")
    		save(algo,algo.opts["filename"])
    	else
    		@warn "could not find `filename` in algo.opts - not saving"
    	end

    	@info "Done with estimation after $t1 minutes"

    	if get(algo.opts,"animate",false)
    		gif(algo.anim,joinpath(dirname(@__FILE__),"../../proposals.gif"),fps=2)
    	end
    end

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
end
