

# MOpt.jl: Moment Optimization Library for Julia

This package implements several MCMC algorithms to optimize a non-differentiable objective function, which depends on a set of parameter values and a set of data moments. In general, this implements *Simulated Method of Moments*. The library is targeted to run MCMC on an SGE cluster, where each node is a chain.

For an R implementation which is the basis of this package, see [https://github.com/tlamadon/mopt](https://github.com/tlamadon/mopt)

## API

### Objective Function

you need an objective function that has an input/output structure as follows:

* Input: parameter vector `p` 
	1. this could be a simple `Array{Float}`. However we want to keep parameter names throughout.
	2. `p` could be a custom type. we access the members of `p` with 
	
	```julia
	i = "alpha"
	newval = 3.4
	eval(parse("p.$(i) = newval"))
	```

* Output: 
	1. function value (model deviation from data moments - metric up to the user): `Float64`
	2. custom model info (user specific. could be a `Dict`)

### Initiating the `Mopt` Object

set up a `MOpt` object by calling the constructor. a list of compulsory arguments could be

* `params_to_sample`: a `DataFrame` with as many rows as parameters you want to sample. There **must** be three columns called `name` (name of param to sample), `lb` (lower bound) and `ub` (upper bound)
* `moments`: a `DataFrame` with 3 columns, `name`, `data` and `sd`. This is the complete set of moments that you have in the data / that is produced by your model
* `moments_to_use`: a vector of moment names you want to use in the current estimatiom. This subsets `moments`.
* your model must produce all moments listed in `moments[:name]`
* your objective function looks something like 

```julia
distance(model(p),moments[findin(moments[:name],moments_to_use),:data])
```

where `distance` is your metric. In other words, your objective is to compute a distance between model and data for the moments in `moments_to_use`.
* `objfunc`: the name of your objective function. Ideally you should export that function so it can be used by `Mopt`
* `mode`: `"serial"` or `"mpi"`


### Running the estimation

call `runMOpt(mymopt)`.






