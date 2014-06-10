
module Mopt


# imports
using  DataFrames
import Base.show


# includes
include("mprob.jl")
include("chains.jl")
include("incmopt.jl")
include("AlgoAbstract.jl")




	
# exports
export MProb, Chain, Testobj
	






end 	# module




