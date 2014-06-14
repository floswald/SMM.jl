
module Mopt


# imports
using  DataFrames, Distributions, Debug
import Base.show


# includes
include("mprob.jl")
include("chains.jl")
include("incmopt.jl")
include("AlgoAbstract.jl")
include("AlgoBGP.jl")




	
# exports
export MProb, Chain, Testobj, getindex
	






end 	# module




