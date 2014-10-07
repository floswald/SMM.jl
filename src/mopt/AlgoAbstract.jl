




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

	# TODO
	# currently this recreates the full dataset on each save
	# if this is slow or otherwise problematic, do something like

    # ff5 = h5open(filename, "w")
	# g = ff5["parameters"]
	# a_data = d_create(g, "a", datatype(Float64), dataspace(algo["maxiter"],1), "chunk", (ifloor(algo["maxiter"] / algo["save_frequency"]),1))
	# b_data = d_create(g, "b", datatype(Float64), dataspace(algo["maxiter"],1), "chunk", (ifloor(algo["maxiter"] / algo["save_frequency"]),1))
	# if saving...
	# a_data[i-algo["save_frequency"] + 1 : i] = ...
	# although needs to be done by chain...

	# do iteration
	for i in 1:algo["maxiter"]

		t0 = time()

		algo.i = i

		computeNextIteration!( algo )

		# save at certain frequency
		# TODO look into HDF5 chunk saving?
		if haskey(algo.opts,"save_frequency")
			@assert haskey(algo.opts,"filename")
			if mod(i,algo["save_frequency"]) == 0
				save(algo,algo["filename"])
			end
		end

		# printing progress
		t1 = round(time()-t0,2)
		if get(algo.opts,"print_level",0) > 2 
			info("iteration $i took $t1 seconds")
			println()
			println(infos(algo.MChains,algo.i))
		elseif get(algo.opts,"print_level",0) > 1
			if mod(algo.i,10) == 0
				info("iteration $i took $t1 seconds")
				println()
				println(infos(algo.MChains,algo.i))
			end
		elseif get(algo.opts,"print_level",0) > 0
			if mod(algo.i,100) == 0
				info("iteration $i took $t1 seconds")
				println()
				println(infos(algo.MChains,algo.i))
			end
		end
	end
	info("Done with estimation loop.")
end

function ps2s_names(algo::MAlgo)
	return ps2s_names(algo.m)
end

function ms_names(algo::MAlgo)
	return ms_names(algo.m)
end

function parameters(m::MAlgo, ch :: Int64, iter:: Int64, p::Symbol)
	return m.MChains[ch].parameters[iter,p]
end

function moments(m::MAlgo, ch :: Int64, iter:: Int64, p::Symbol)
	return m.MChains[ch].moments[iter,p]
end

function moments(m::MAlgo, ch :: Int64, iter:: Int64)
	return [ m.MChains[ch].moments[iter,p] for p in ms_names(m)]
end

