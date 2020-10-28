

module SMM

__precompile__(false)

# Dependencies
# ############

using Distributions
using DataFrames, DataFramesMeta
using OrderedCollections
using PDMats
using Documenter
using Plots
using FileIO
using JLD2
using StatsPlots
# using JSON
using ProgressMeter
using Random
using Distributed
using LinearAlgebra
using Logging
using Statistics



import Base.get, Base.write, Base.==
import Base.show, Statistics.mean, Statistics.median

# exports: Types
export MProb, Eval, MAlgo, MAlgoBGP 

# exports: methods
export addParam!,
       addSampledParam!,
       addMoment!,
       addEvalFunc!,
       setMoments!,
       setValue!,
       readEvalArray ,
       dataMoment,
       dataMomentd,
       dataMomentW,
       summary,
       params,
       param,
       paramd,
       doSlices,
       start,
       finish,
       save,
       readMalgo,
       restart!,
       extendBGPChain!,
       run!,
       mean,median,CI

# create a random device that is immune to seeds
const RAND = RandomDevice()


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
include("mopt/econometrics.jl")


end 	# module
