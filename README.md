

# MOpt.jl: Moment Optimization Library for [Julia](http://julialang.org)


[![Build Status](https://travis-ci.org/floswald/MOpt.jl.png?branch=master)](https://travis-ci.org/floswald/MOpt.jl)

This package provides a `Julia` infrastructure for *[Simulated Method of Moments](http://en.wikipedia.org/wiki/Method_of_simulated_moments)* estimation, or other problems where we want to optimize a non-differentiable objective function. The setup is suitable for all kinds of **likelihood-free estimators** - in general, those require evaluating the objective at many regions. The user can supply their own algorithms for generating successive new parameter guesses. We provide a set of MCMC template algorithms. The code can be run in serial or on a cluster.

[![acceptance rates](doc/img/acceptance.png)]()

## Detailed Documentation

[Documentation available on readthedocs.](http://moptjl.readthedocs.org/en/latest/)


## Example Usage of the BGP Algorithm

```julia
using MOpt

# data are generated from a bivariate normal
# with mu = [a,b] = [0,0], and cov = Diagonal([1,1])
# aim: 
# 1) sample [a',b'] from a space [-2,2] x [-1,1] and
# 2) find true [a,b] by computing distance(S([a',b']), S([a,b]))
#    and accepting/rejecting [a',b'] according to BGP
# 3) S([a,b]) returns a summary of features of the data

# initial value
p    = ["a" => 1.9 , "b" => -0.9]
# param bounds
pb   = [ "a" => [-2,2] , "b" => [-1,1] ]
# data moments
moms = DataFrame(moment=["alpha","beta"],data_value=[0.0,0.0],data_sd=rand(2))

# define a minization problem
mprob = MProb(p,pb,MOpt.objfunc_norm2,moms)

# look at slices of the model: 
# how do objective function and
# simulated moments vary with each
# parameter one by one, holding 
# the others fixed at initial value?

obj_slices = MOpt.slices(mprob,30)
MOpt.plotSlices(mprob,obj_slices[1],obj_slices[2])
```

[![objective slices](doc/img/slices_objective.png)]()
[![alpha slices](doc/img/slices_alpha.png)]()
[![beta slices](doc/img/slices_beta.png)]()


```julia
# setup a minization algorithm: options

opts =[
	"N"               => 6,							# number of MCMC chains
	"mode"            => "serial",					# mode: serial or mpi
	"maxiter"         => 500,						# max number of iterations
	"savefile"        => joinpath(pwd(),"MA.h5"),	# filename to save results
	"print_level"     => 1,							# increasing verbosity level of output
	"maxtemp"         => 100,						# tempering of hottest chain
	"min_shock_sd"    => 0.1,						# initial sd of shock on coldest chain
	"max_shock_sd"    => 1,							# initial sd of shock on hottest chain
	"past_iterations" => 30,						# num of periods used to compute Cov(p)
	"min_accept_tol"  => 100,						# ABC-MCMC cutoff for rejecting small improvements
	"max_accept_tol"  => 100,						# ABC-MCMC cutoff for rejecting small improvements
	"min_disttol"     => 0.1,						# distance tol for jumps from coldest chain
	"max_disttol"     => 0.1,						# distance tol for jumps from hottest chain
	"min_jump_prob"   => 0.05,						# prob of jumps from coldest chain
	"max_jump_prob"   => 0.2]						# prob of jumps from hottest chain

# setup the BGP algorithm
MA = MAlgoBGP(mprob,opts)

# run it
runMopt!(MA)

# plot outputs
plot(MA,"acc")
# see first plot above
```

```julia
plot(MA,"params_time")
```

[![acceptance rates](doc/img/pars_time.png)]()

```julia
plot(MA,"params_dist")
```

[![acceptance rates](doc/img/pars_dist.png)]()

```julia
# save results
save(MA,MA["savefile"])
```

## Contributing

We encourage user contributions. Please submit a pull request for any improvements you would like to suggest, or a new algorithm you implemented. 

New algorithms:
*You can model your algo on the basis of `src/AlgoBGP.jl` - 
* your algorithm needs a storage system for chains, derived from the default type `Chain` in `src/chains.jl`
* you need to implement the function `computeNextIteration!( algo )` for your `algo`










