using Random
using LinearAlgebra
using Statistics
using Distributions
using DataFrames
using OrderedCollections
using SMM

Random.seed!(1234)

function datagenOLS(N=10_000,K=4)
    X = 10*randn(N,K)
    X[:,1] .= 1
    β = [2, -1, 1, 0.5]
    σ = .5
    θ = vcat(β,σ)
    ε = σ*randn(N)
    y = X*β .+ ε
    return y,X,θ
end

function objfunc_ols(ev::Eval,Xd)
    start(ev)
    
    # in OLS: dataMoment is mean( y.*X[:,k] ) for all k
    # in OLS: simMoment  is mean( (X*beta .+ eps).*X[:,k] ) for all k

    # extract parameters
    θ  = collect(values(ev.params))
    nm = length(ev.dataMoments)

    # compute simulated moments
    if get(ev.options,:noseed,false)
    else
        Random.seed!(1234)
    end
    ns = 10_000
    N = size(Xd,1)
    K = size(Xd,2)
    β = θ[1:end-1]
    σ = θ[end]
    ε = sqrt(σ^2)*randn(N) # take sqrt of σ^2 to handle case where optimizer guesses a negative value for σ
    simM = zeros(K+1)
    for k = 1:K
        simM[k] = mean( (Xd*β .+ ε).*Xd[:,k] )
    end
    simM[K+1] = var(Xd*β) + σ^2 # this is var(y) in the simulation, and we want it to match var(y) in the data

    # get data moments
    # same thing here, and use dataMomentsWeights for sd
    # second argument can be optional
    # get objective value: (data[i] - model[i]) / weight[i]
    v = Dict{Symbol,Float64}()
    simMoments = Dict{Symbol,Float64}()
    i = 0
    for (k,mom) in dataMomentd(ev)
        i += 1
        simMoments[k] = simM[i]
        if haskey(dataMomentWd(ev),k)
            v[k] = ((simMoments[k] .- mom) ./ dataMomentW(ev,k)) .^2
        else
            v[k] = ((simMoments[k] .- mom) ) .^2
        end
    end
    setValue!(ev, mean(collect(values(v))))

    # also return the moments
    setMoments!(ev, simMoments)

    ev.status = 1

    # finish and return
    finish(ev)

    return ev
end

function parallelOLS(startvals,svlb,svub,y,X,niter=200)
    # data are generated from a bivariate normal
    # with μ = [a,b] = [0,0]
    # aim:
    # 1) sample [a',b'] from a space [-3,3] x [-2,2] and
    # 2) find true [a,b] by computing distance(S([a',b']), S([a,b]))
    #    and accepting/rejecting [a',b'] according to BGP
    # 3) S([a,b]) returns a summary of features of the data: 2 means

    K = size(startvals,1)-1 # number of X's in model (including intercept)
    # starting values
    pb = OrderedDict()
    for k = 0:K-1
        pb[string("b",k)] = [startvals[k+1],svlb[k+1],svub[k+1]]
    end
    pb["s"] = [startvals[K+1],svlb[K+1],svub[K+1]]
    # parameter names
    names = String[]
    for k = 0:K-1
        push!(names,string("β",k))
    end
    push!(names,"σ") # SD of epsilon
    # parameter estimates
    values = Float64[]
    for k = 1:K
        push!(values,mean( y.*X[:,k] ))
    end
    push!(values,var(y)) # want to match var(y) in simulation and actual data to be able to estimate σ

    # pass these objects to MProb
    moms = DataFrame(name=names,value=values,weight=ones(K+1))
    println(moms)
    mprob = MProb()
    addSampledParam!(mprob,pb)
    addMoment!(mprob,moms)
    addEvalFunc!(mprob,objfunc_ols) # how to pass X into `objfunc_ols`??

    nchains = 3

    opts =Dict("N"=>nchains,
        "maxiter"=>niter,
        "maxtemp"=> 5,
            # choose inital sd for each parameter p
            # such that Pr( x \in [init-b,init+b]) = 0.975
            # where b = (p[:ub]-p[:lb])*opts["coverage"] i.e. the fraction of the search interval you want to search around the initial value
        "coverage"=>0.025,  # i.e. this gives you a 95% CI about the current parameter on chain number 1.
        "sigma_update_steps"=>10,
        "sigma_adjust_by"=>0.01,
        "smpl_iters"=>100, # 1000
        "parallel"=>false,
        "min_improve"=>[0.05 for i in 1:nchains],
        "mixprob"=>0.3,
        "acc_tuners"=>[12.0 for i in 1:nchains],
        "animate"=>false)

    # setup the BGP algorithm
    MA = MAlgoBGP(mprob,opts)

    # run the estimation
    run!(MA)
    @show SMM.summary(MA)

    return MA
end

y,X,theta = datagenOLS()
println(X\y)
sval_lb = theta .- round.(2 .* abs.(theta);digits=1)
sval_ub = theta .+ round.(2 .* abs.(theta);digits=1)
MA = parallelOLS(rand(length(theta)),sval_lb,sval_ub,y,X)
dc = SMM.history(MA.chains[1])
dc = dc[dc[:accepted].==true, :]
println(describe(dc))
