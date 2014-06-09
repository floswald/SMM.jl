
# there are 2 types of information that will be required. 
#
# 1. information about the prbolem to solve, such as
# the number of parameters, the moments, etc ...  this will be stored in MProb
#
# 2. information about the way to solve the problem, such as which algorithm to use,
# the parameters for this algorithm, the paths, etc... this will be stored in MAlgo





# tried that: doesnt work because the Mopt needs to take a generic parameter type, like a vector or a dict. otherwise every user will have their own param type. I don't even know the name of that user type, let alone it's structure.
# a test param type
# type Testparam
# 	a :: Float64
# 	b :: Float64

# 	function Testparam(a,b) 
# 	return new(a,b)
# 	end
# end



# call this function to set up
# the cluster. 
# maps choice of "mode" into an action
function MoptPrepare!(m::Moptim)

	if m.prepared
		println("you're good to go")
		return nothing
	else

		# setup cluster
		@everywhere include(m.include_on_workers)

		# set m.N = # of workers

		# send data to workers

		# set m.prepared = true
	end

end


function evaluateObjective(m::Moptim)

	if m.mode =="serial"
		v = map(x -> evaluateChainID(x,m), 1:m.N )
	else
		v = pmap(x -> evaluateChainID(x,m), 1:m.N )

	end
	return v

end

function evaluateChainID(m::Moptim,i::Int)

	# eval chain i with param[i]
	x = eval(Expr(:call,m.objfunc,m.current_param[i],m.moments,m.moments_to_use))
	return x

end

function collectChains(m::Moptim)



end






# I think was very nice below.
# evalChain! operates inplace, we don't have to
# go and collect results from an array

# # evaluate the objective function
# # on all chains at their respective
# # parameter values
# function evaluateObjective(m::Moptim)

# 	if m.mode=="serial"
# 		for v in values(m.chains) 
# 			evalChain!(v,m) 
# 		end
# 	else
# 		if !m.prepared
# 			error("you must call MoptPrepare before evaluating the objective")
# 		end
# 		pmap( x -> evalChain!(x,m) , values(m.chains))
# 	end
# 	return nothing
# end




# # evalute chain
# # calls the objective function with the 
# # parameter vector currently stored in that chain
# function evalChain!(c::MCMChain,m::Moptim)

# 	# call objective function of m
# 	# x could be a tuple of values: (value, time, status)
# 	x = eval(Expr(:call,m.objfunc,c.p,m.moments,m.moments_to_use))

# 	# TODO check length of x
# 	# and hope the user returns the correct order?

# 	# store results in c
# 	# 
# 	# make a vector out of dict p
# 	tmp = DataFrame(id=c.id,iter=c.iter,value=x[1],time=x[2],status=x[3])
# 	# add param
# 	tmp = hcat(tmp,)
# 	c.data = rbind(c.data,)
# 	c.phist[c.iter] = c.p

# 	return nothing
# end

# update param on a single chain
# function updateChain!(c::MCMChain,newp::Dict)
# 	c.p = newp
# 	c.iter += 1
# 	return nothing
# end

# # update all chains
# function updateAllChains!(m::Moptim,newps::Dict)
# 	for x in keys(m.chains) 
# 		updateChain!(m.chains[x],newps[x]) 
# 	end
# 	return nothing
# end

# # update param
# function updateParam(m::Moptim)

# 	# for all chains in m.chains
# 	# compute a new guess

# 	# this is the actual MCMC algorithm
# 	# for each chain takes current and last value
# 	# and accepts/rejects it

# end



