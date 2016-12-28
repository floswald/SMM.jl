
module MOpt

# Dependencies
# ############

using Distributions 
using Logging
using DataFrames, DataFramesMeta
using HDF5
import Base.show, Base.std
using GLM
using DataStructures
using PDMats
using Roots

@Logging.configure(level=DEBUG)


# exports: Types
export MProb, Eval,MAlgo, MAlgoBGP 
# export MProb, Chain, BGPChain, Testobj, Eval, AbstractChain

# exports: methods
export addParam!, addSampledParam!, addMoment!, addEvalFunc!, setMoment, start, finish, param, paramd , setValue, readEval,write, readEvalArray , dataMoment,dataMomentW
#        dataMomentd, dataMomentWd,setMoment
# export getindex, setindex, parameters, evals, infos, allstats, moments, hist,
       # runMOpt!, save, slices, transpose, , readEvalArrayRemote, 


# load files
# ############

include("mopt/mprob.jl")
include("mopt/incmopt.jl")
include("mopt/Eval.jl")
# include("mopt/chains.jl")
# include("mopt/slices.jl")
include("mopt/AlgoAbstract.jl")
include("mopt/AlgoBGP.jl")
include("mopt/ObjExamples.jl")
# include("mopt/Examples.jl")
# include("mopt/sobolsearch.jl")
# include("mopt/econometrics.jl")

# if is_apple()
#        using PyPlot
#        include("mopt/plotting.jl")
# end


end 	# module




