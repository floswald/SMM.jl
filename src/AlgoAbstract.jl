




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




function runMopt!( algo::MAlgo )

	# tasks
	# =====
	# load data from file if set in algo.opts
	# setup cluster if required

	# do iteration
	for i in 1:algo["maxiter"]

		algo.i = i

		computeNextIteration!( algo )

		# save at certain frequency
		# TODO look into HDF5 chunk saving

		if mod(algo.i,100) == 0
			println("algo.i=$(algo.i)")
		end
	end
end








function histogram(x,nb::Int) 
  n, bins = hist(x,nb)
  nvec = [n]
  bar(nvec[1:end-1], bins, width = nvec[2] - nvec[1])
end

function plot(algo::MAlgo, what)

	if (what == "acc")
		dd = infos(algo.MChains)
		subplot(211)
		for sdf in groupby(dd, :chain_id)
  			PyPlot.plot(sdf[:iter],sdf[:accept_rate])
  			I = sdf[:exchanged_with].>0
  			PyPlot.plot( sdf[I,:iter], sdf[I,:accept_rate],"o")
		end
		subplot(212)
		for sdf in groupby(dd, :chain_id)
  			PyPlot.plot(sdf[:shock_sd])
  			I = sdf[:exchanged_with].>0
  			PyPlot.plot( sdf[I,:iter], sdf[I,:shock_sd],"o")
		end
		suptitle("Acceptance Rates by Chain")
	end

	if (what == "params_time")
		dd = sort!(parameters(algo.MChains,true)[[:chain_id,:iter,algo.m.p2sample_sym]])
		npars = length(algo.m.p2sample_sym)
		nrows = floor(sqrt(npars))
		ncols = ceil(npars/nrows)
		pid = 0
		for par in algo.m.p2sample_sym
			pid += 1
			subplot(nrows,ncols,pid)
			for sdf in groupby(dd, :chain_id)
  				plot(sdf[:iter],sdf[par])
  			end
		end
		suptitle("Parameter values over time")
	end
	if (what == "params_dist")
		dd = parameters(algo.MChains,true)[[:chain_id,:iter,algo.m.p2sample_sym]]
		npars = length(algo.m.p2sample_sym)
		nrows = floor(sqrt(npars))
		ncols = ceil(npars/nrows)
		pid = 0
		for par in algo.m.p2sample_sym
			pid += 1
			subplot(nrows,ncols,pid)
			# println(sdf[:iter])
			# plot(sdf[:iter],sdf[par])
			histogram(dd[par],30)
			title(string(par))
			xlabel("parameter value")
  		end
		suptitle("Posterior Distribution of Parameters")
	end

end





