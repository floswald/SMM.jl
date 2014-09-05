Getting Started
===============

Let's see how to get started very quickly. You start by ....


Setting the problem
-------------------

This package implements several MCMC algorithms to optimize a non-differentiable objective function. The main application are **likelihood-free estimators**, which requires evaluating the objective at many regions. In general, this implements *Simulated Method of Moments*. The library is targeted to run MCMC on an SGE cluster, where each node is a chain.

For an R implementation which is the basis of this package, see [https://github.com/tlamadon/mopt](https://github.com/tlamadon/mopt)

Example Useage
--------------

.. code-block:: julia

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


  
