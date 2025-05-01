



"""
This abstract type nests all [`MProb`](@ref) algorithms, for example [`SMM.MAlgoBGP`](@ref)
"""
abstract type MAlgo end

import Base.getindex, Base.setindex!

# getter and setters for Algo.opts
function getindex(algo::MAlgo, key)
  return(algo.opts[key])
end

function setindex!(algo::MAlgo, val,key)
  algo.opts[key] = val
end


"""
	run!( algo::MAlgo )

Function to start estimation of an [`MAlgo`](@ref).
"""
function run!( algo::MAlgo )

	# tasks
	# =====
	# load data from file if set in algo.opts
	# setup cluster if required

	@info "Starting estimation loop."
	t0 = time()

	# do iteration
    @showprogress for i in 1:algo["maxiter"]
    # @showprogress 1 "Running Estimation..." for i in 1:algo["maxiter"]
		@debug "iteration $i"

		algo.i = i

		# try
			computeNextIteration!( algo )

			# save at certain frequency
			if haskey(algo.opts,"save_frequency") == true
				if haskey(algo.opts,"filename") == true
	  				if mod(i,algo.opts["save_frequency"]) == 0
	  					save(algo,algo.opts["filename"])
	  					# @info "saved data at iteration $i")
	  				end
        		end
			end

		# catch e
		# 	@warn "caught exception $e")
		# 	throw(e)
		# end
	end
	t1 = round((time()-t0)/60,digits = 1)
	algo.opts["time"] = t1
	if haskey(algo.opts,"filename")
		save(algo,algo.opts["filename"])
	else
		@warn "could not find 'filename' in algo.opts - not saving anything"
	end

	@info "Done with estimation after $t1 minutes"

	if get(algo.opts,"animate",false)
		gif(algo.anim,joinpath(dirname(@__FILE__),"../../proposals.gif"),fps=2)
	end

end

"""
	save(algo::MAlgo, filename::AbstractString)

Save MAlgo to disk using JLD2
"""
function save(algo::MAlgo, filename::AbstractString)

   # saving the entire MAlgoBGP object:
   JLD2.@save filename algo

end

"""
	readMalgo(filename::AbstractString)

Load MAlgo from disk
"""
function readMalgo(filename::AbstractString)

    # load the entire MAlgoBGP object:
    JLD2.@load filename algo

    return algo

end

