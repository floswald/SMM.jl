
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
    best_id   :: Vector{Int}   # index of best eval.value so far
    best_val  :: Vector{Float64}   # best eval.value so far
    curr_val  :: Vector{Float64}   # current value
    probs_acc :: Vector{Float64}    # vector of probabilities with which to accept
    mprob     :: MProb
    id        :: Int64
    iter      :: Int64
    accepted  :: Array{Bool}
    accept_rate :: Float64
    acc_tuner :: Float64
    exchanged :: Array{Int}
    m         :: MProb
    sigma     :: PDiagMat{Float64}
    sigma_update_steps :: Int64   # update sampling vars every sigma_update_steps iterations
    sigma_adjust_by :: Float64   # adjust sampling vars by sigma_adjust_by percent up or down
    smpl_iters :: Int64   # max number of trials to get a new parameter from MvNormal that lies within support
    maxdist  :: Float64  # what's the maximal function value you will accept when proposed a swap. i.e. if ev.value > maxdist, you don't want to swap with ev.

    function BGPChain(id::Int=1,n::Int=10,m::MProb=MProb(),sig::Vector{Float64}=Float64[],upd::Int64=10,upd_by::Float64=0.01,smpl_iters::Int=1000,maxdist::Float64=10.0,acc_tuner::Float64=2.0)
        @assert length(sig) == length(m.params_to_sample)
        this           = new()
        this.mprob     = m
        this.evals     = Array{Eval}(n)
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
        this.sigma     = PDiagMat(sig)
        this.sigma_update_steps = upd
        this.sigma_adjust_by = upd_by
        this.smpl_iters = smpl_iters
        this.maxdist = maxdist
        return this
    end
end

allAccepted(c::BGPChain) = c.evals[c.accepted]

# return a dict of param values as arrays
function params(c::BGPChain)
    e = allAccepted(c)
    d = Dict{Symbol,Vector{Float64}}()
    for k in keys(e[1].params)
        d[k] = Float64[e[i].params[k] for i in 1:length(e)]
    end
    return d
end

# make a dataframe with acc/rej history
function history(c::BGPChain)
    N = length(c.evals)
    cols = Any[]
    # d = DataFrame([Int64,Float64,Bool,Int64],[:iter,:value,:accepted,:prob],N)
    d = DataFrame()
    d[:iter] = collect(1:c.iter)
    d[:exchanged] = c.exchanged
    d[:accepted] = c.accepted
    d[:best_val] = c.best_val
    d[:curr_val] = c.curr_val
    d[:best_id] = c.best_id
    # get fields from evals
    nms = [:value,:prob]
    for n in nms
        d[n] = eltype(getfield(c.evals[1],n))[getfield(c.evals[i],n) for i in 1:N]
    end
    # get fields from evals.params
    for (k,v) in c.evals[1].params
        d[k] = eltype(v)[c.evals[i].params[k] for i in 1:N]
    end

    return d[[:iter,:value,:accepted,:curr_val, :best_val, :prob, :exchanged,collect(keys(c.evals[1].params))...]]
end

best(c::BGPChain) = findmin([c.evals[i].value for i in 1:length(c.evals)])
mean(c::BGPChain) = Dict(k => mean(v) for (k,v) in params(c))

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
        return find(c.accepted[1:(c.iter)])[end]
    end
end
getIterEval(c::BGPChain,i::Int) = c.evals[i]
getLastAccepted(c::BGPChain) = c.evals[lastAccepted(c)]
set_sigma!(c::BGPChain,s::Vector{Float64}) = length(s) == length(c.m.params_to_sample) ? c.sigma = PDiagMat(s) : ArgumentError("s has wrong length")
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
            c.best_id[c.iter]  = c.iter
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

function next_eval!(c::BGPChain)
    # generate new parameter vector from last accepted param

    # increment interation
    c.iter += 1
    @debug("iteration = $(c.iter)")

    # returns an OrderedDict
    pp = proposal(c)

    # evaluate objective 
    ev = evaluateObjective(c.m,pp)

    # accept reject 
    doAcceptReject!(c,ev)

    # save eval on BGPChain 
    set_eval!(c,ev)

end
# function Remote_next_eval!(conn,pp::OrderedDict)
#     # generate new parameter vector from last accepted param

#     write(conn, pp)

#     # other stuff needs to be moved to master
#     # ev = evaluateObjective(c.m,pp)

#     # # accept reject 
#     # doAcceptReject!(c,ev)

#     # # save eval on BGPChain 
#     # set_eval!(c,ev)

# end



function doAcceptReject!(c::BGPChain,eval_new::Eval)
            @debug("")
            @debug("doAcceptReject!")
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

            # this forumulation: old - new
            # because we are MINIMIZING the value of the objective function            
            eval_new.prob = minimum([1.0,exp( c.acc_tuner * ( eval_old.value - eval_new.value) )]) #* (eval_new.value < )
            @debug("eval_new.value = $(eval_new.value)")
            @debug("eval_old.value = $(eval_old.value)")
            @debug("eval_new.prob = $(round(eval_new.prob,2))")
            @debug("c.probs_acc[c.iter] = $(round(c.probs_acc[c.iter],2))")

            if isna(eval_new.prob) || !isfinite(eval_new.prob)
                eval_new.prob = 0.0
                eval_new.accepted = false
                eval_new.status = -1

            elseif !isfinite(eval_old.value)
                # should never have gotten accepted
                @debug("eval_old is not finite")
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
            @debug("eval_new.accepted = $(eval_new.accepted)")
            @debug("")
            @debug("")

        end

        c.accepted[c.iter] = eval_new.accepted
        set_acceptRate!(c)

        # update sampling variances every x periods
        # -----------------------------------------

        # update shock variance. want to achieve a long run accpetance rate of 23.4% (See Casella and Berger)
        # and only if you are not BGPChain number 1

        # if (c.id>1) && (mod(c.iter,c.sigma_update_steps) == 0)

        # if mod(c.iter,c.sigma_update_steps) == 0
        #     too_high = c.accept_rate > 0.234
        #     if too_high
        #         @debug("acceptance rate on BGPChain $(c.id) is too high at $(c.accept_rate). increasing variance of each param by $(100* c.sigma_adjust_by)%.")
        #         set_sigma!(c,diag(c.sigma) .* (1.0+c.sigma_adjust_by) )
        #     else
        #         @debug("acceptance rate on BGPChain $(c.id) is too low at $(c.accept_rate). decreasing variance of each param by $(100* c.sigma_adjust_by)%.")
        #         set_sigma!(c,diag(c.sigma) .* (1.0-c.sigma_adjust_by) )
        #     end
        # end
    end
end


"""
    mysample(d::Distributions.MultivariateDistribution,lb::Vector{Float64},ub::Vector{Float64},iters::Int)

mysample from distribution `d` until all poins are in support
"""
function mysample(d::Distributions.MultivariateDistribution,lb::Vector{Float64},ub::Vector{Float64},iters::Int)

    # draw until all points are in support
    for i in 1:iters
        x = rand(d)
        if all(x.>=lb) && all(x.<=ub)
            return x
        end
    end
    error("no draw in support after $iters trials. increase either opts[smpl_iters] or opts[bound_prob].")
end

function proposal(c::BGPChain)

    if c.iter==1
        return c.m.initial_value
    else
        ev_old = getLastAccepted(c)
        mu  = paramd(ev_old) # dict of params
        lb = [v[:lb] for (k,v) in c.m.params_to_sample]
        ub = [v[:ub] for (k,v) in c.m.params_to_sample]

        # Transition Kernel is q(.|theta(t-1)) ~ TruncatedN(theta(t-1), Sigma,lb,ub)
        newp = Dict(zip(collect(keys(mu)),mysample(MvNormal(collect(values(mu)),c.sigma),lb,ub,c.smpl_iters)))
        # @debug("iteration $(c.iter)")
        # @debug("old param: $(ev_old.params)")
        # @debug("new param: $newp")

        # flat kernel: random choice in each dimension.
        # newp = Dict(zip(collect(keys(mu)),rand(length(lb)) .* (ub .- lb))) 

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
  
    function MAlgoBGP(m::MProb,opts=Dict("N"=>3,"maxiter"=>100,"maxtemp"=> 2,"coverage"=>0.125,"sigma_update_steps"=>10,"sigma_adjust_by"=>0.01,"smpl_iters"=>1000,"parallel"=>false,"maxdists"=>[0.5 for i in 1:3],"mixprob"=>0.5,"acc_tuner"=>2))

        init_sd = OrderedDict{Symbol,Float64}()
        if opts["N"] > 1
    		temps     = linspace(1.0,opts["maxtemp"],opts["N"])
            # initial std dev for each parameter to achieve at least bound_prob on the bounds
            # println("opts=$opts")
            # println("pars = $( m.params_to_sample)")

            # choose inital sd for each parameter p
            # such that Pr( x \in [init-b,init+b]) = 0.975
            # where b = (p[:ub]-p[:lb])*opts["coverage"] i.e. the fraction of the search interval you want to search around the initial value
            for (k,v) in m.params_to_sample
                # mu = (v[:lb]+v[:ub])/2
                b = (v[:ub]-v[:lb])*opts["coverage"]
                # init_sd[k] = MOpt.initsd(mu+b,mu)
                # @assert init_sd[k] == b / quantile(Normal(),0.975)
                init_sd[k] = b / quantile(Normal(),0.975)
            end
            BGPChains = BGPChain[BGPChain(i,opts["maxiter"],m,collect(values(init_sd)) .* temps[i],opts["sigma_update_steps"],opts["sigma_adjust_by"],opts["smpl_iters"],opts["maxdists"][i],opts["acc_tuner"]) for i in 1:opts["N"]]
          else
            temps     = [1.0]
            for (k,v) in m.params_to_sample
                b = (v[:ub]-v[:lb])*opts["coverage"]
                init_sd[k] = b / quantile(Normal(),0.975)
            end
            # println(init_sd)
            BGPChains = BGPChain[BGPChain(1,opts["maxiter"],m,collect(values(init_sd)),opts["sigma_update_steps"],opts["sigma_adjust_by"],opts["smpl_iters"],opts["maxdists"][1],opts["acc_tuner"])]
        end
	    return new(m,opts,0,BGPChains)
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

#     r = Dict()
#     for ic in 1:length(m.chains)
#         if m.i == 0
#             r[ic] = Dict(:mu => m.m.initial_value,:sigma => m.chains[ic].sigma)
#         else
#             ev_old = getLastAccepted(m.chains[ic])
#             r[ic] = Dict(:mu => paramd(ev_old),:sigma => m.chains[ic].sigma)
#         end
#     end
#     r
# end

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




# computes new candidate vectors for each BGPChain
# accepts/rejects that vector on each BGPChain, according to some rule
# *) computes N new parameter vectors
# *) applies a criterion to accept/reject any new params
# *) stores the result in BGPChains
function computeNextIteration!( algo::MAlgoBGP )
    # here is the meat of your algorithm:
    # how to go from p(t) to p(t+1) ?

    # incrementBGPChainIter!(algo.chains)


    # TODO 
    # this on each BGPChain
    # START=========================================================
    if get(algo.opts, "parallel", false) 
        pmap( x->next_eval!(x), algo.chains ) # this does proposal, evaluateObjective, doAcceptRecject
    else
        # for i in algo.chains
        #     @debug(" ")
        #     @debug(" ")
        #     @debug("debugging chain id $(i.id)")
        #     next_eval!(i)
        # end
        map( x->next_eval!(x), algo.chains ) # this does proposal, evaluateObjective, doAcceptRecject
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


function exchangeMoves!(algo::MAlgoBGP)

    # algo["N"] exchange moves are proposed
    props = [(i,j) for i in 1:algo["N"], j in 1:algo["N"] if (j<i)]
    # N pairs of chains are chosen uniformly in all possibel pairs with replacement
    samples = algo["N"] < 3 ? algo["N"]-1 : algo["N"]
    pairs = sample(props,samples,replace=false)

    @debug("")
    @debug("exchangeMoves: proposing pairs")
    @debug("$pairs")

    for p in pairs
        i,j = p
        evi = getLastAccepted(algo.chains[p[1]])
        evj = getLastAccepted(algo.chains[p[2]])
        if rand() < algo["mixprob"]
            # exchange i with j if rho(S(z_j),S(data)) < epsilon_i
            if (evj.value < evi.value)  # if j's value is better than i's
            # if (abs(j.value) < algo.chains[p[1]].best_val) && (rand() < algo["mixprob"]) # notice: tolerance on chain i!
                @debug("$j better than $i")
                # @debug("$(abs(j.value)) < $(algo.chains[p[1]].maxdist)")
                # swap_ev!(algo,p)
                set_ev_i2j!(algo,i,j)
            else
                @debug("$i better than $j")
                set_ev_i2j!(algo,j,i)
            end
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
 #                    @debug("perc dist $ch and $ch2 is $tmp. will label that `close`.")
	# 				push!(close,ch2)
	# 			end
	# 		end
	# 	end
	# 	# 2) with y% probability exchange with a randomly chosen BGPChain from close
	# 	if length(close) > 0
	# 		ex_with = rand(close)
	# 		@debug("making an exchange move for BGPChain $ch with BGPChain $ex_with set: $close")
	# 		swap_ev!(algo,Pair(ch,ex_with))
	# 	end
	# end

end

function swap_ev!(algo::MAlgoBGP,ch_id::Tuple{Int,Int})
    @debug("swapping evs for $(ch_id[1]) and $(ch_id[2])")
    c1 = algo.chains[ch_id[1]] 
    c2 = algo.chains[ch_id[2]]

	e1 = getLastAccepted(c1)
	e2 = getLastAccepted(c2)

    # swap 
    set_eval!(c1,e2)
    set_eval!(c2,e1)

	# make a note 
    set_exchanged!(c1,ch_id[2])
    set_exchanged!(c2,ch_id[1])
end
function set_ev_i2j!(algo::MAlgoBGP,i::Int,j::Int)
    @debug("setting ev of $i to ev of $j")
    ci = algo.chains[i] 
    cj = algo.chains[j]

    ei = getLastAccepted(ci)
    ej = getLastAccepted(cj)

    # set ei -> ej
    set_eval!(ci,ej)

    # make a note 
    set_exchanged!(ci,j)
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
end
