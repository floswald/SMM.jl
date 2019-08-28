




# A Moment Optimizing Algorithm is called MAlgo
abstract type MAlgo end

import Base.getindex, Base.setindex!

# getter and setters for Algo.opts
function getindex(algo::MAlgo, key)
  return(algo.opts[key])
end

function setindex!(algo::MAlgo, val,key)
  algo.opts[key] = val
end

function runMOpt!( algo::MAlgo )

	# tasks
	# =====
	# load data from file if set in algo.opts
	# setup cluster if required

	@info(logger,"Starting estimation loop.")
	t0 = time()

	# do iteration
	println("running estimation")
    for i in 1:algo["maxiter"]
    	print("*")
    # @showprogress 1 "Running Estimation..." for i in 1:algo["maxiter"]
		@debug(logger,"iteration $i")

		algo.i = i

		# try
			computeNextIteration!( algo )

			# save at certain frequency
			if haskey(algo.opts,"save_frequency") == true
				if haskey(algo.opts,"filename") == true
	  				if mod(i,algo.opts["save_frequency"]) == 0
	  					save(algo,algo.opts["filename"])
	  					# @info(logger,"saved data at iteration $i")
	  				end
        		end
			end

		# catch e
		# 	@warn(logger,"caught exception $e")
		# 	throw(e)
		# end
	end
	t1 = round((time()-t0)/60,1)
	algo.opts["time"] = t1
	if haskey(algo.opts,"filename")
		save(algo,algo.opts["filename"])
	else
		@warn(logger,"could not find 'filename' in algo.opts")
	end

	@info(logger,"\nDone with estimation after $t1 minutes")

	if get(algo.opts,"animate",false)
		gif(algo.anim,joinpath(dirname(@__FILE__),"../../proposals.gif"),fps=2)
	end

end

"""
  save(algo::MAlgoBGP, filename::AbstractString)

Save MAlgo to disk using JLD2
"""
function save(algo::MAlgo, filename::AbstractString)

   # saving the entire MAlgoBGP object:
   JLD2.@save filename algo

end

"""
   load(filename::AbstractString)

Load MAlgo from disk
"""
function readMalgo(filename::AbstractString)

    # load the entire MAlgoBGP object:
    JLD2.@load filename algo

    return algo

end



# function ps2s_names(algo::MAlgo)
# 	return ps2s_names(algo.m)
# end

# function ms_names(algo::MAlgo)
# 	return ms_names(algo.m)
# end

# function parameters(m::MAlgo, ch :: Int64, iter:: Int64, p::Symbol)
# 	return m.MChains[ch].parameters[iter,p]
# end

# function moments(m::MAlgo, ch :: Int64, iter:: Int64, p::Symbol)
# 	return m.MChains[ch].moments[iter,p]
# end

# function moments(m::MAlgo, ch :: Int64, iter:: Int64)
# 	return [ m.MChains[ch].moments[iter,p] for p in ms_names(m)]
# end
