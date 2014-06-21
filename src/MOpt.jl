
module MOpt

# Dependencies
# ############

using Distributions, PyPlot
using Reexport
@reexport using DataFrames
import Base.show
import PyPlot.plot

# exports
# ############

export MProb, 
       Chain,  
       BGPChain,  
       MAlgo,
       MAlgoBGP,
       Testobj, 
       getindex,
       setindex,
       parameters,
       evals,
       infos,
       allstats,
       moments,
       hist,
       runMopt!


# load files
# ############

include("mprob.jl")
include("chains.jl")
include("incmopt.jl")
include("AlgoAbstract.jl")
include("AlgoBGP.jl")


end 	# module




