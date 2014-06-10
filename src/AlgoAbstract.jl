



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
	algo.MChains.chains[which]
end

function updateIter!(algo::MAlgo)
	algo.i =+ 1
	return nothing
end


# updates all chains in an algo
# wrapper to evaluateChainID
# calls map(evaluateChainID) or pmap(evaluateChainID) depending on opts
function updateChains!(algo::MAlgo)

	# add 1 to iteration count on algo
	updateIter!(algo)
	computeNextIteration!( algo )

end
    

# evalute objective function
# with param vector number i
function evaluateObjective(algo::MAlgo,which::Int)

	# eval chain i with param p
	x = eval(Expr(:call,algo.m.objfunc,algo.candidate_param[which],algo.m.moments,algo.m.moments_subset))
	return x

end




