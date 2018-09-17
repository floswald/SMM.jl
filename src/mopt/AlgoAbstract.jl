




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

	@info "Starting estimation loop."
	t0 = time()

	# do iteration
	for i in 1:algo["maxiter"]
		@debug(logger,"iteration $i")

		algo.i = i

		# try
			computeNextIteration!( algo )

			# save at certain frequency
			if haskey(algo.opts,"save_frequency") == true
        # if the user provided a filename in the options dictionary
				if haskey(algo.opts,"filename") == true
  				if mod(i,algo.opts["save_frequency"]) == 0
  					save(algo,algo.opts["filename"])
  					@info "saved data at iteration $i"
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
    # if no filename is provided, generated a random number
    	filename,err = mktemp()
		@warn "could not find 'filename' in algo.opts"
    	save(algo,filename)
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
   tempfilename = filename * ".jld2"
   JLD2.@save tempfilename algo

end

"""
   load(filename::AbstractString)

Load MAlgo from disk
"""
function readMalgo(filename::AbstractString)

    # load the entire MAlgoBGP object:
    tempfilename =  filename * ".jld2"
    JLD2.@load tempfilename algo

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
