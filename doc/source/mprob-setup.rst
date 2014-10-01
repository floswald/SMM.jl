Setting up a moment problems
=============================

This section describes how to create the `MProb` object that will store the description of the problem. 
In general this is composed of the set of parameters we will iterate over together with their bounds, as well as the set of moments to match, with their value and their precision.

We start by creating an empty MProb object and we progressively add content to it.

.. code-block:: julia

  using Mopt
  mprob = MProb()

We then describe the step required to scuccessfuly set up the object.

Step 1: add parameters
----------------------

There are two types of parameters, some are sampled by the algorythm and some are not. Here are ways to add arguments to the description:

.. code-block:: julia

  addParam!(mprob, "c", 0.1)
  addParam!(mprob, "d", 0.2)
  addSampledParam!(mprob, "a", 0.1, 0, 1)
  addSampledParam!(mprob, "b", 0.1, 0, 1)

But parameters can also be added all at once using for instance a dictionary

.. code-block:: julia

  ps  = { "c" => 0.1 , "d" => 0.2}
  addParam!(mprob,ps)
  pss = { "a" => [0.1, 0, 1] , "b" => [0.1, 0, 1]}
  addSampledParam!(mprob,pss)

Step 2: add moments
-------------------

The second step is to add moments to the description. These are moments that the objective function will use as matching criteria. This is more for analysis purposes, as ex-post it can be very informative to look at what parameter is affected by what particular moments, or to check the shape of the moment conditions at maximum.

In a spirit similar to the parameters, moments can be added as follows:

.. code-block:: julia

  addMoment!(mprob, "m1", 0.1, 0.001)
  addMoment!(mprob, "m2", 0.1, 0.001)

But as well using a DataFrame which is usually how moments are loaded, from a csv file or other source. In this case the used needs to provide the names of the columns that must be used.

.. code-block:: julia

  dd = DataFrame(name= ["m1", "m2"], value=[0.1,0.1], sd=[0.01,0.01])
  addMoment!(mprob, dd, [:name,:value,:sd])

Step 3: set the objective function
----------------------------------

This is simply giving the function that will evaluate the objective for a given set of parameters. For this example, we use a test objective function defined within the package. See the dedicated section on supplying your own obejective function.

.. code-block:: julia

  addObj!(mprob, objfunc_norm2 )

Done, next is selecting an algorithm
------------------------------------

This is it, your problem is now created, you can move to the second step which is to select an algorithm.

Putting all at once using Lazy.jl
---------------------------------

.. code-block:: julia

  mprob = @> begin
    Mprob()
    addParam!(       "c", 0.1)
    addParam!(       "d", 0.2)
    addSampledParam!("a", 0.1, 0, 1)
    addSampledParam!("b", 0.1, 0, 1) 
    addMoment!( "m1", 0.1, 0.001)
    addMoment!( "m2", 0.1, 0.001)
    addObj!( objfunc_norm2 )
  end


