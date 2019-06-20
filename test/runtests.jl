

using MomentOpt
using Test
using DataFrames
using OrderedCollections
using Distributed
using LinearAlgebra
using Random
using Plots: Animation
using JSON

# If we want the test to pass, we need this
# see https://github.com/JuliaPlots/Plots.jl/issues/1076
# Otherwise the following error shows up:
# "GKS: Open failed in routine OPEN_WS
# GKS: GKS not in proper state. GKS must be either in the state WSOP or WSAC in routine ACTIVATE_WS
# /home/travis/.julia/v0.6/GR/src/../deps/gr/bin/gksqt: error while loading shared libraries: libQt5Widgets.so.5: cannot open shared object file: No such file or directory
# connect: Connection refused
# GKS: can't connect to GKS socket application
# Did you start 'gksqt'?""
ENV["GKSwstype"] = "100"

include(joinpath(dirname(@__FILE__),"include","test-include.jl"))

include("test_MProb.jl")
include("test_Eval.jl")
include("test_slices.jl")
include("test_BGPChain.jl")
include("test_AlgoAbstract.jl")
include("test_AlgoBGP.jl")
include("test_objfunc.jl")


