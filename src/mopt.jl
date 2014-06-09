
module Mopt


# imports
using  DataFrames
import Base.show


# includes
include("chains.jl")
include("incmopt.jl")

# MCMChain methods
function summaryChain(c::MCMChain)

	# compute summary stats of chain

end



# exports
export Moptim, Testparam
	






end 	# module




