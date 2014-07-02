
require("src/nodes.jl")


# compute slices in parallel
mprob_slices = MOpt.slices(mprob,30)
println(mprob_slices[1])
println(mprob_slices[2])



