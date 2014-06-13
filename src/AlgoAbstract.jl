



# define an abstract type and set/get for it
abstract MAlgo

# getter and setters for Algo.opts
function getindex(algo::MAlgo, key)
  return(algo.opts[key])
end

function setindex!(algo::MAlgo, val,key)
  algo.opts[key] = val
end

# get chain number "which" from algo
function getChain(algo::MAlgo, which::Int)
	algo.chains[which]
end

function updateIter!(algo::MAlgo)
	algo.i =+ 1
	return nothing
end


# updates all chains in an algo
# wrapper to evaluateChainID
# calls map(evaluateChainID) or pmap(evaluateChainID) depending on opts
function updateChains!(algo::MAlgo)
	
	updateIter!(algo)  # add 1 to iteration count on algo
	computeNextIteration!( algo )  # compute next iteration on all chains

end
    

# evalute objective function
# with param vector number i
# always evaluates the field "candidate_param"
function evaluateObjective(algo::MAlgo,which::Int)

	# eval chain i with param p
	x = eval(Expr(:call,algo.m.objfunc,algo.candidate_param[which],algo.m.moments,algo.m.moments_subset))
	return x

end




