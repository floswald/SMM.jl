
module MomentOpt

__precompile__(false)

# Dependencies
# ############

using Logging
using Distributions 
using DataFrames, DataFramesMeta
import Base.show 
using Statistics
using DataStructures
using PDMats
using Documenter
using Plots
using FileIO
using JLD2
using GLM
using Distributed

gr()


# set log level
LogLevel(Logging.Debug)

import Base.get,Base.write, Base.start, Base.==
import Statistics: mean

# exports: Types
export MProb, Eval,MAlgo, MAlgoBGP 

# exports: methods
export addParam!, 
       addSampledParam!, 
       addMoment!, 
       addEvalFunc!, 
       setMoment, 
       setValue,
       readEvalArray , 
       dataMoment,
       dataMomentd,
       dataMomentW,
       summary,
       param,
       paramd,
       doSlices,
       start,
       finish,
       save,
       readMalgo,
       restartMOpt!,
       extendBGPChain!,
       fitMirror



# load files
# ############

include("mopt/mprob.jl")
include("mopt/incmopt.jl")
include("mopt/Eval.jl")
include("mopt/slices.jl")
include("mopt/AlgoAbstract.jl")
include("mopt/AlgoBGP.jl")
include("mopt/ObjExamples.jl")
include("mopt/Examples.jl")
include("mopt/plotting.jl")


end    # module




