




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




function runMOpt!( algo::MAlgo )

	# tasks
	# =====
	# load data from file if set in algo.opts
	# setup cluster if required

	info("Starting estimation loop.")

	# do iteration
	for i in 1:algo["maxiter"]

		algo.i = i

		computeNextIteration!( algo )

		# save at certain frequency
		# TODO look into HDF5 chunk saving

		if haskey(algo.opts,"print_level")
			if algo["print_level"] > 2
				println(infos(algo.MChains,algo.i))
			elseif algo["print_level"] > 1
				if mod(algo.i,10) == 0
					println(infos(algo.MChains,algo.i))
				end
			elseif algo["print_level"] > 0
				if mod(algo.i,100) == 0
					println(infos(algo.MChains,algo.i))
				end
			end
		end
	end
end







