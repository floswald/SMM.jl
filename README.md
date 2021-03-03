

# SMM.jl: Simulated Method of Moments for Julia

| **Documentation** | **Build Status**                                                                                |
|:---------------------:|:------------------:|
| [![][docs-stable-img]][docs-stable-url] [![][docs-dev-img]][docs-dev-url] | [![][travis-img]][travis-url] [![][appveyor-img]][appveyor-url] |

*Notice*: this package was previously called *MomentOpt.jl*.


This package provides a `Julia` infrastructure for *[Simulated Method of Moments](http://en.wikipedia.org/wiki/Method_of_simulated_moments)* estimation, or other problems where we want to optimize a non-differentiable objective function. The setup is suitable for all kinds of **likelihood-free estimators** - in general, those require evaluating the objective at many regions. The user can supply their own algorithms for generating successive new parameter guesses. We provide a set of MCMC template algorithms. The code can be run in serial or on a cluster.


## Installation

In your julia REPL, type

```julia
] add SMM
```

## Documentation

[![][docs-stable-img]][docs-stable-url]


### Examples

Please check out a fully worked example in [`src/Examples.jl`](src/Examples.jl).

Here is a working session comparing serial vs parallel performance on a test objective function. Notice that parallel performance hinges on the objective function being reasonably expensive to compute (at least 0.1 seconds per function evaluation) - otherwise the overhead from data transfer is just too high.

```julia
julia> using SMM
[ Info: Precompiling SMM [bc769cb7-f249-5bba-802a-79f18cb247ec]

julia> x = SMM.serialNormal(2,200,slow = true)
[ Info: Starting estimation loop.
Progress: 100%|██████████████████████████████████████████████████████████████████████████████████████████████| Time: 0:01:05
┌ Warning: could not find 'filename' in algo.opts - not saving anything
└ @ SMM ~/.julia/dev/SMM/src/mopt/AlgoAbstract.jl:67
[ Info: Done with estimation after 1.1 minutes
summary(MA) = 3×5 DataFrame
 Row │ id     acc_rate   perc_exchanged  exchanged_most_with  best_val
     │ Int64  Float64    Float64         Int64                Float64
─────┼──────────────────────────────────────────────────────────────────
   1 │     1  0.0793651             5.5                    2  0.0023224
   2 │     2  0.0819672             8.5                    1  0.0126754
   3 │     3  0.115183              4.5                    2  0.0145625
(
BGP Algorithm with 3 BGPChains
============================

Algorithm
---------
Current iteration: 200
Number of params to estimate: 2
Number of moments to match: 2

, Plot{Plots.GRBackend() n=2}, Plot{Plots.GRBackend() n=3})

julia> using Distributed

julia> addprocs(2, exeflags="--project=.")  # you don't need the `exeflag` if you `add`ed the package regularly!
2-element Array{Int64,1}:
 2
 3

julia> @everywhere using SMM

julia> x = SMM.serialNormal(2,200,slow = true)
[ Info: Starting estimation loop.
Progress: 100%|██████████████████████████████████████████████████████████████████████████████████████████████| Time: 0:00:49
┌ Warning: could not find 'filename' in algo.opts - not saving anything
└ @ SMM ~/.julia/dev/SMM/src/mopt/AlgoAbstract.jl:67
[ Info: Done with estimation after 0.8 minutes
summary(MA) = 3×5 DataFrame
 Row │ id     acc_rate   perc_exchanged  exchanged_most_with  best_val
     │ Int64  Float64    Float64         Int64                Float64
─────┼───────────────────────────────────────────────────────────────────
   1 │     1  0.117347              2.0                    2  0.00246371
   2 │     2  0.0899471             5.5                    3  0.103399
   3 │     3  0.161458              4.0                    2  0.139263
(
BGP Algorithm with 3 BGPChains
============================

Algorithm
---------
Current iteration: 200
Number of params to estimate: 2
Number of moments to match: 2

, Plot{Plots.GRBackend() n=2}, Plot{Plots.GRBackend() n=3})
```

## Contributing

We encourage user contributions. Please submit a pull request for any improvements you would like to suggest, or a new algorithm you implemented.

New algorithms:
* You can model your algo on the basis of `src/AlgoBGP.jl` -
* you need to implement the function `computeNextIteration!( algo )` for your `algo`

## History

This package grew out of the [R package `mopt`](https://github.com/tlamadon/mopt). 

## Thanks to all Contributors!

* [Julien Pascal](https://github.com/JulienPascal)
* [Edoardo Ciscato](https://github.com/edoardociscato)

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://floswald.github.io/SMM.jl/latest

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://floswald.github.io/SMM.jl/latest

[travis-img]: https://travis-ci.org/floswald/SMM.jl.svg?branch=master
[travis-url]: https://travis-ci.org/floswald/SMM.jl

[appveyor-img]: https://ci.appveyor.com/api/projects/status/github/floswald/SMM.jl?branch=master&svg=true
[appveyor-url]: https://ci.appveyor.com/project/floswald/SMM.jl/branch/master
