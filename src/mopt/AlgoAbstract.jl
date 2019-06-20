



"""
# Moment Minisation Algorithm Type

This abstract type nests all [`MProb`](@ref) algorithms, for example [`AlgoBGP`](@ref)
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
	runMOpt!( algo::MAlgo )

Function to start estimation of an [`MAlgo`]@(ref).
"""
function runMOpt!( algo::MAlgo )

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
		@warn "could not find 'filename' in algo.opts"
	end

	@info "Done with estimation after $t1 minutes"

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
