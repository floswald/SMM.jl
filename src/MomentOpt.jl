
module MomentOpt

# Dependencies
# ############

using Distributions 
using MiniLogging
using DataFrames, DataFramesMeta
using HDF5
import Base.show, Base.std
using DataStructures
using PDMats
using Documenter
using Plots

pyplot

# setup MiniLogging
logger = get_logger()
if isinteractive()
    basic_config(MiniLogging.DEBUG; date_format="%H:%M:%S")
else
    basic_config(MiniLogging.INFO; date_format="%H:%M:%S")
end

import Base.get, Base.mean, Base.write, Base.start

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
       finish



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


end 	# module




