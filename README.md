

# MOpt.jl: Moment Optimization Library for Julia

This package implements several MCMC algorithms to optimize a non-differentiable objective function. The main application are **likelihood-free estimators**, which requires evaluating the objective at many regions. In general, this implements *Simulated Method of Moments*. The library is targeted to run MCMC on an SGE cluster, where each node is a chain.

For an R implementation which is the basis of this package, see [https://github.com/tlamadon/mopt](https://github.com/tlamadon/mopt)

## Example Useage

```julia
using Mopt

# get a parameter vector
p = ["a" => 3.1 , "b" => 4.9]
# define params to use with bounds
pb= [ "a" => [0,1] , "b" => [0,1] ]

# get some moments
# first entry is moment estimate, second is standard deviation
moms = [
	"alpha" => [ 0.8 , 0.02 ],
	"beta"  => [ 0.8 , 0.02 ],
	"gamma" => [ 0.8 , 0.02 ]
]

# a subset of moments to match
submoms = ["alpha", "beta"]

# call objective
x = Mopt.Testobj(p,moms,submoms)

# Define an Moment Optimization Problem
mprob = Mopt.MProb(p,pb,Mopt.Testobj,moms;moments_subset=submoms)

# show
mprob

# step 2: choose an algorithm
# ----------------------
algo = Mopt.MAlgoRandom(mprob,opts=["mode"=>"serial","maxiter"=>100])

# step 3: run estimation
# ----------------------
runMopt(algo)

```

## API

* there are three custom types: `MProb`, `MAlgo` and `MChain`

### `MProb`

#### Objective Function

You supply an objective function to the constructor of `MProb` must be callable as `objfunc(p,moments,whichmom,...)`, where the first three arguments are compulsory and are explained below.

you need an objective function that has an input/output structure as follows:

Input: 

1. parameter vector `p`: parameterized as a `Dict`. Thought about custom types (hard to implement because I'd have to guess what exactly that type looks like), and simple (unnamed) vectors. Took `Dict` because closest to `R list()`
2. data moments. 
3. `whichmom`: `Array{ASCIIString}` of moment names to which you want the objective function subset to in the current estimation

Output: 
1. function value (model deviation from data moments - metric up to the user): `Float64`
2. custom model info (user specific. could be a `Dict`)

#### Initiating the `MProb` Object

set up a `MProb` object by calling the constructor. there is a list of **compulsory arguments**, which must be supplied in the following **order**: 

1. `p`: the initial value of the full parameter vector.
2. `params_to_sample`: a `Dict{ASCIIString,Array{Float,1}}` of the kind `["a" => [0,1]]`, where `a` is the name of the parameter and `[0,1]]` are lower and upper bounds. (all  parameters not in `params_to_sample` remain fixed). 
3. `objfunc`: type `Function` representing your objective function. Internally we evaluate `Expr(:call,m.objfunc,m.current_param,m.moments,m.moments_to_use)`. 
4. `moments`: a `Dict{ASCIIString,Array{Float,1}}` ["alpha" => [1.192,0.01]]`, where `alpha` is the name of the moment and `[1.192,0.01]]` are value and standard deviation.

Then there is a list of **optional arguments** which have to be supplied by *keyword*. If you don't supply them, we'll take a default value:

* `moments_to_use`: a vector of moment names you want to use in the current estimatiom. This subsets `moments`. Default=all


### MAlgo

this type implements most of the computation. There is a lot of flexibility to define new algorithms. 

options:

* `shock_var`: variance of shocks. Default=1 
* `np_shock`: new parameter shock. Default=1
* `save_freq`: after how many iterations to save to disk. Default=25
* `N`: Number of chains. If `mode=="mpi"`, choose `N==numberOfNodes`. Default=3
* `n_untempered`: number of untempered chains. Default =1
* `mode`: "serial" or "mpi". Default="serial"
* `path`: a string to where to save  







