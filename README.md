

# MOpt.jl: Moment Optimization Library for [Julia](http://julialang.org)

[![Build Status](https://travis-ci.org/floswald/MOpt.jl.png?branch=master)](https://travis-ci.org/floswald/MOpt.jl)

This package provides a `Julia` infrastructure for *[Simulated Method of Moments](http://en.wikipedia.org/wiki/Method_of_simulated_moments)* estimation, or other problems where we want to optimize a non-differentiable objective function. The setup is suitable for all kinds of **likelihood-free estimators** - in general, those require evaluating the objective at many regions. The user can supply their own algorithms for generating successive new parameter guesses. We provide a set of MCMC template algorithms. The code can be run in serial or on a cluster.

**Track proposals by iteration**  

![Poposals](https://github.com/floswald/MOpt.jl/master/proposals.gif)

## Installation

```julia
Pkg.clone("https://github.com/floswald/MOpt.jl")
```

## Documentation

Still work in progress. :-(  
I recommend to look at `src/Examples.jl`:

### Example Usage of the BGP Algorithm

Baragatti, Grimaud and Pommeret (BGP) in ["Likelihood-free parallel tempring"](http://arxiv.org/abs/1108.3423) propose an approximate Bayesian Computation (ABC) algorithm that incorporates the parallel tempering idea of Geyer (1991). We provide the BGP algorithm as a template called `MAlgoBGP`. Here we use it to run a simple toy example.

```julia
using MOpt

# this is MOpt.serialNormal()

# data are generated from a bivariate normal
# with mu = [a,b] = [0,0]
# aim: 
# 1) sample [a',b'] from a space [-3,3] x [-2,2] and
# 2) find true [a,b] by computing distance(S([a',b']), S([a,b]))
#    and accepting/rejecting [a',b'] according to BGP
# 3) S([a,b]) returns a summary of features of the data: 2 means

# initial value
pb    = OrderedDict("p1" => [0.2,-3,3] , "p2" => [-0.2,-2,2] )
moms = DataFrame(name=["mu1","mu2"],value=[-1.0,1.0],weight=ones(2))
mprob = MProb() 
addSampledParam!(mprob,pb) 
addMoment!(mprob,moms) 
addEvalFunc!(mprob,objfunc_norm)

nchains = 3

opts =Dict("N"=>nchains,
	"maxiter"=>200,
	"maxtemp"=> 5,
        # choose inital sd for each parameter p
        # such that Pr( x \in [init-b,init+b]) = 0.975
        # where b = (p[:ub]-p[:lb])*opts["coverage"] i.e. the fraction of the search interval you want to search around the initial value
	"coverage"=>0.025,  # i.e. this gives you a 95% CI about the current parameter on chain number 1.
	"sigma_update_steps"=>10,
	"sigma_adjust_by"=>0.01,
	"smpl_iters"=>1000,
	"parallel"=>false,
	"maxdists"=>[0.05 for i in 1:nchains],
	"mixprob"=>0.3,
	"acc_tuner"=>12.0,
	"animate"=>true)

# setup the BGP algorithm
MA = MAlgoBGP(mprob,opts)

# run the estimation
runMOpt!(MA)
@show summary(MA)

│ Row │ id │ acc_rate │ perc_exchanged │ exchanged_most_with │ best_val   │
├─────┼────┼──────────┼────────────────┼─────────────────────┼────────────┤
│ 1   │ 1  │ 0.581967 │ 39.0           │ 2                   │ 1.85755e-5 │
│ 2   │ 2  │ 0.346939 │ 27.0           │ 3                   │ 1.85755e-5 │
│ 3   │ 3  │ 0.245399 │ 18.5           │ 2                   │ 1.85755e-5 │
```

### Example Output

** param histogram **  

![Histogram](https://github.com/floswald/MOpt.jl/master/histogram.png)  


** param history **  

![Lines](https://github.com/floswald/MOpt.jl/master/lines.png)  

** Slices of objective function **  

![Slices](https://github.com/floswald/MOpt.jl/master/slices-v.png)  

** param histogram **  

![Slices1](https://github.com/floswald/MOpt.jl/master/slices-m.png)  
![Slices2](https://github.com/floswald/MOpt.jl/master/slices-m2.png)  


## Contributing

We encourage user contributions. Please submit a pull request for any improvements you would like to suggest, or a new algorithm you implemented. 

New algorithms:
*You can model your algo on the basis of `src/AlgoBGP.jl` - 
* your algorithm needs a storage system for chains, derived from the default type `Chain` in `src/chains.jl`
* you need to implement the function `computeNextIteration!( algo )` for your `algo`










