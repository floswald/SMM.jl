# SMM.jl

*Simulated Method of Moments (SMM) for Julia*

A package providing supporting infrastructure and algorithms to perform [Simulated Method of Moments](http://en.wikipedia.org/wiki/Method_of_simulated_moments).

## Features

* Fully fledged SMM infrastructure to track moments and function evaluations.
* Supply your own objective function to be optimized.
* Optimization can be carried out in parallel.
* [Diagnostic tools](@ref) illustrating the proposal distribution as well as markov chain statistics.
* Includes a parellel tempering [likelihood-free optimizer](http://arxiv.org/abs/1108.3423) using MCMC technology and a [coordinate descent optimizer](https://en.wikipedia.org/wiki/Coordinate_descent).

For some example usage see the [Examples](@ref) page.


## Manual Outline

```@contents
Pages = ["index.md","eval.md","slices.md","algo.md","metrics.md","examples.md"]
```


## `MProb`: Minimisation/Maximisation Problems

* A core object of this library is the `MProb` type, specifying an optimisation problem.

```@autodocs
Modules = [SMM]
Pages = ["mprob.jl"]
```
