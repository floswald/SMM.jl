




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

function runMopt!( algo::MAlgo )

	# tasks
	# =====

	# load data from file if set in algo.opts

	# setup cluster if required

	# do iteration
	for i in 1:algo["maxiter"]

		algo.i = i

		if mod(algo.i,100) == 0
			println("algo.i=$(algo.i)")
		end

		computeNextIteration!( algo )

		# save at certain frequency

		# reporting


	end

	# save
end

function hist(x,nb=30) 
  n, bins = PyPlot.hist(x, nb)
  bar(n[1:end-1], bins, width = n[2] - n[1])
end

function plot(algo::MAlgo, what)

	if (what == "acc")
		dd = Mopt.infos(algo.MChains)
		subplot(211)
		for sdf in Mopt.groupby(dd, :chain_id)
  			PyPlot.plot(sdf[:iter],sdf[:accept_rate])
  			I = sdf[:exchanged_with].>0
  			PyPlot.plot( sdf[I,:iter], sdf[I,:accept_rate],"o")
		end
		subplot(212)
		for sdf in Mopt.groupby(dd, :chain_id)
  			PyPlot.plot(sdf[:shock_sd])
  			I = sdf[:exchanged_with].>0
  			PyPlot.plot( sdf[I,:iter], sdf[I,:shock_sd],"o")
		end
	end

	if (what == "params_time")
		dd = sort!(Mopt.parameters(algo.MChains))
		npars = length(algo.m.p2sample_sym)
		nrows = floor(sqrt(npars))
		ncols = ceil(npars/nrows)
		pid = 0
		# for par in algo.m.p2sample_sym
		# 	pid += 1
		# 	subplot(nrows,ncols,pid)
		# 	for ic in 1:algo["N"]
		# 		# println(sdf[:iter])
  # 				# plot(sdf[:iter],sdf[par])
  # 				plot(algo.MChains[ic].parameters[:iter],algo.MChains[ic].parameters[par])
  # 			end
		# end
		for par in algo.m.p2sample_sym
			pid += 1
			subplot(nrows,ncols,pid)
			for sdf in Mopt.groupby(dd, :chain_id)
  				plot(sdf[:iter],sdf[par])
  			end
		end
	end
	if (what == "params_dist")
		dd = Mopt.parameters(algo.MChains)
		npars = length(algo.m.p2sample_sym)
		nrows = floor(sqrt(npars))
		ncols = ceil(npars/nrows)
		pid = 0
		for par in algo.m.p2sample_sym
			pid += 1
			subplot(nrows,ncols,pid)
			# println(sdf[:iter])
			# plot(sdf[:iter],sdf[par])
			hist(dd[par])
			title(string(par))
			xlabel("parameter value")
  		end
	end

end





