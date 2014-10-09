
module MOpt

# Dependencies
# ############

using Distributions, HDF5
using Reexport
using Lumberjack
@reexport using DataFrames
import Base.show, Base.transpose

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
       runMOpt!,
       save,
       slices,
       transpose,
       addParam!,
       addSampledParam!


# load files
# ############

include("mopt/mprob.jl")
include("mopt/incmopt.jl")
include("mopt/Eval.jl")
include("mopt/chains.jl")
include("mopt/slices.jl")
include("mopt/AlgoAbstract.jl")
include("mopt/AlgoBGP.jl")
include("mopt/ObjExamples.jl")
include("mopt/sobolsearch.jl")


# for now plotting only on my box because
# installing matplotlib on unix hpc is tricky.
# comment this out if you want the plot function
if Sys.OS_NAME == :Darwin
       using PyPlot
       include("mopt/plotting.jl")
end


end 	# module




