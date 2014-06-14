



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



# evalute objective function
# with param vector number i
# always evaluates the field "candidate_param"
function evaluateObjective(algo::MAlgo,which::Int)

	# eval chain i with param p
	x = eval(Expr(:call,algo.m.objfunc,algo.candidate_param[which],algo.m.moments,algo.m.moments_subset))
	return x

end

function runMopt( algo::MAlgo )

	# tasks
	# =====

	# load data from file if set in algo.opts

	# setup cluster if required

	# do iteration
	for i in 1:algo["maxiter"]

		algo.i = i

		computeNextIteration!( algo )

		# save at certain frequency

		# reporting

	end

	# save
end



