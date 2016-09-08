
module MOpt

# Dependencies
# ############

using Distributions 
using DataFrames
using HDF5
import Base.show, Base.std
using GLM


# exports: Types
export MProb, Chain, BGPChain, MAlgo, MAlgoBGP, Testobj, Eval, AbstractChain

# exports: methods
export getindex, setindex, parameters, evals, infos, allstats, moments, hist,
       runMOpt!, save, slices, transpose, addParam!, addSampledParam!, addMoment!, addEvalFunc!, start, finish, param, paramd, 
       dataMoment,dataMomentd, dataMomentW, dataMomentWd,setMoment, setValue, readEval, readEvalArray, readEvalArrayRemote, write


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
include("mopt/Examples.jl")
include("mopt/sobolsearch.jl")
include("mopt/econometrics.jl")


# for now plotting only on my box because
# installing matplotlib on unix hpc is tricky.
# comment this out if you want the plot function
if Sys.OS_NAME == :Darwin
       using PyPlot
       include("mopt/plotting.jl")
end

end 	# module




