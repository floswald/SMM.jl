
module MOpt

# Dependencies
# ############

using Distributions 
using Reexport
using Lumberjack
@reexport using DataFrames
import Base.show
if VERSION < v"0.4.0-dev"
    using Docile
end
@document

"This package implements several MCMC algorithms to optimize a non-differentiable objective function. The main application are **likelihood-free estimators**, which requires evaluating the objective at many regions. In general, this implements *Simulated Method of Moments*. The library is targeted to run MCMC on an SGE cluster, where each node is a chain."

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
       addSampledParam!,
       addMoment!,
       Eval,
       start,
       finish,
       param,
       paramd,
       fill,
       dataMoment,
       dataMomentW,
       setMoment,
       setValue,
       readEval,
       readEvalArray,
       readEvalArrayRemote,
       write

if !haskey(ENV,"IGNORE_HDF5")
       using HDF5
end

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




