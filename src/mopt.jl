
module Mopt


# imports
using  DataFrames
import Base.show


# includes
include("mprob.jl")
include("chains.jl")
include("incmopt.jl")
include("AlgoAbstract.jl")
include("Algos.jl")




	
# exports
export MProb, Chain, Testobj, getindex
	






end 	# module




