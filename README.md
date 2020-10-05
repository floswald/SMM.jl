# Forked version of SMM.jl

This is a modified version of SMM.jl with two changes for the STBNews project:

1. Change to parallel dispatch to improve memory usage (in `algoBGP.jl`). Note STB code uses `algoSTB.jl`, in `stb-model` project.
2. Allow passing options to objective function (in `mprob.jl::evaluateObjective`)

To install, in your julia REPL, type

```julia
] add https://github.com/gregobad/SMM.jl#better_parallel
```


# SMM.jl: Simulated Method of Moments for Julia

| **Documentation** | **Build Status**                                                                                |
|:---------------------:|:------------------:|
| [![][docs-stable-img]][docs-stable-url] [![][docs-dev-img]][docs-dev-url] | [![][travis-img]][travis-url] [![][appveyor-img]][appveyor-url] |

*Notice*: this package was previously called *MomentOpt.jl*.


This package provides a `Julia` infrastructure for *[Simulated Method of Moments](http://en.wikipedia.org/wiki/Method_of_simulated_moments)* estimation, or other problems where we want to optimize a non-differentiable objective function. The setup is suitable for all kinds of **likelihood-free estimators** - in general, those require evaluating the objective at many regions. The user can supply their own algorithms for generating successive new parameter guesses. We provide a set of MCMC template algorithms. The code can be run in serial or on a cluster.

## Documentation

[![][docs-stable-img]][docs-stable-url]


### Example Notebook

Please check out a fully worked example in [`src/example.ipynb`](src/example.ipynb).

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
