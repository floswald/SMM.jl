
module MOpt

# Dependencies
# ############

using Distributions, PyPlot
using Reexport
@reexport using DataFrames
import Base.show, Base.transpose
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
       savePlots,
       transpose


# load files
# ############

include("mopt/mprob.jl")
include("mopt/chains.jl")
include("mopt/incmopt.jl")
include("mopt/AlgoAbstract.jl")
include("mopt/AlgoBGP.jl")


end 	# module




