

using MomentOpt
using Base.Test
using TestSetExtensions
using DataFrames

include(joinpath(dirname(@__FILE__),"include","test-include.jl"))

@testset ExtendedTestSet "Running MomentOpt tests" begin

    @includetests ARGS

end
    



