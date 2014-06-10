



# define an abstract type and set/get for it
abstract MAlgo

# getter and setters for Algo.opts
function getindex(algo::MAlgo, key)
  return(algo.opts[key])
end

function setindex!(algo::MAlgo, val,key)
  algo.opts[key] = val
end

function computeNextIteration!( algo::MAlgo  )
  error("computeNextIteration not implemented for this algorithm")
end

# get chain number "which" from algo
function getChain(algo::MAlgo, which::Int)
	algo.chains.chains[which]
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

	if algo.i>1 
		# get a new parameter vector for each chain
		computeNextIteration!( algo )
	end

end
    

# evalute chain number i
# with param vector number i
function evaluateChainID(algo::MAlgo,i::Int)

	# eval chain i with param[i]
	x = eval(Expr(:call,algo.m.objfunc,algo.current_param[i],algo.m.moments,algo.m.moments_subset))
	return x

end




