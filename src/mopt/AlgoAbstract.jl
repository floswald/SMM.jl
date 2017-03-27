




# A Moment Optimizing Algorithm is called MAlgo
abstract MAlgo

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

	@info("Starting estimation loop.")
	t0 = time()

	# do iteration
	for i in 1:algo["maxiter"]

		algo.i = i

		try
			computeNextIteration!( algo )

			# save at certain frequency
			if haskey(algo.opts,"save_frequency")
				@assert haskey(algo.opts,"filename")
				if mod(i,algo["save_frequency"]) == 0
					save(algo,algo["filename"])
					@info("saved data at iteration $i")
				end
			end

		catch e
			@warn("caught exception $e")
			throw(e)
		end
	end
	t1 = round((time()-t0)/60,1)
	algo.opts["time"] = t1
	if haskey(algo.opts,"filename")
		save(algo,algo["filename"])
	else
		@warn("could not find 'filename' and did not save")
	end
	@info("Done with estimation after $t1 minutes")

	if get(algo.opts,"animate",false)
		gif(algo.anim,joinpath(Pkg.dir("MOpt"),"proposals.gif"),fps=2)
	end

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



