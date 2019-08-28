

using MomentOpt
using Base.Test
using TestSetExtensions
using DataFrames
using DataStructures
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

@testset ExtendedTestSet "Running MomentOpt tests" begin

    @includetests ARGS

end
