
module MOpt

# Dependencies
# ############

using Distributions, PyPlot
using Reexport
@reexport using DataFrames
import Base.show
import PyPlot.plot
using HDF5

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
       plot,
       allstats,
       moments,
       hist,
       runMOpt!,
       save,
       slices,
       plotSlices,
       savePlots


# load files
# ############

include("mprob.jl")
include("chains.jl")
include("incmopt.jl")
include("AlgoAbstract.jl")
include("AlgoBGP.jl")


end 	# module




