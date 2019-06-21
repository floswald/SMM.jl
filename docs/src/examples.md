# Examples


## Estimating Means of a bivariate normal


Let's define a function that returns a Normal distribution with a certain location, and let's call this *our model*:

```math
m(p) = \mathcal{N}\left( [p_1,p_2]' , I_2 \right)
```

Our aim is to estimate the location of means from data that is simulated from that law with an MCMC chain. (Of course the sample mean would be a perfectly valid estimator.) The twist here is that we will pretend that we don't have access to the entire dataset ``\{X_i\}_{i=1}^N, X_i = (x_{i1},x_{i2})``, but only a set of summary statistics ``S`` - in our case, we'd have two moments of this data, namely ``\mu_j = \frac{1}{N}x_{ij},j=1,2``. Our objective function is the squared distance between ``\mu``, and what our model produces instead. That is, we give a parameter vector ``p`` to our model ``m``, which in turn produces 2 *simulated* moments denoted ``\hat{\mu}``. Finally, we assess their respective distance.

Again: 

1. Assume *true* moments ``\mu``.
1. repeatedly create data from ``m(p)`` for different ``p`` drawn from a space [-3,3] x [-2,2]. For each dataset, compute ``\hat{\mu_j} = \frac{1}{N}x_{ij},j=1,2``.
1. Compute distance ``\mu,\hat{\mu}`` and decide according to [`BGP algorithm`](@ref) whether to accept or reject current ``p``. 

```julia-repl
julia> using MomentOpt
julia> pb    = OrderedDict("p1" => [0.2,-3,3] , "p2" => [-0.2,-2,2] )
julia> moms  = DataFrame(name=["mu1","mu2"],value=[-1.0,1.0],weight=ones(2))
julia> mprob = MProb() 
julia> addSampledParam!(mprob,pb) 
julia> addMoment!(mprob,moms) 
julia> addEvalFunc!(mprob,objfunc_norm)

julia> nchains = 3

julia> opts =Dict("N"=>nchains,"maxiter"=> 10,"maxtemp"=> 5,"coverage"=>0.025,
 "sigma_update_steps"=>10, "sigma_adjust_by"=>0.01, "smpl_iters"=>1000,
  "parallel"=>true, "min_improve"=>[0.05 for i in 1:nchains], "mixprob"=>0.3, 
  "acc_tuners"=>[12.0 for i in 1:nchains], "animate"=>false)

julia> MA = MAlgoBGP(mprob,opts)

BGP Algorithm with 3 BGPChains
============================

Algorithm
---------
Current iteration: 0
Number of params to estimate: 2
Number of moments to match: 2

julia> runMOpt!(MA)
[ Info: Starting estimation loop.
Progress: 100%|████████████████████████████████████████| Time: 0:00:04
┌ Warning: could not find 'filename' in algo.opts
└ @ MomentOpt ~/.julia/v0.6/MomentOpt/src/mopt/AlgoAbstract.jl:69
[ Info: Done with estimation after 0.1 minutes

```

Full list of options is available at the [`MomentOpt.BGPChain`](@ref) documentation