

using MomentOpt
using Test
using DataFrames
using Distributed
using Random
using GLM
using DataStructures
using Statistics

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

dir = dirname(@__FILE__)
include(joinpath(dir,"include","test-include.jl"))

@testset "Running MomentOpt tests" begin

    include(joinpath(dir,"test_MProb.jl"))
    include(joinpath(dir,"test_Eval.jl"))
    # include(joinpath(dir,"test_slices.jl"))
    # include(joinpath(dir,"test_AlgoAbstract.jl"))
    include(joinpath(dir,"test_AlgoBGP.jl"))

end
