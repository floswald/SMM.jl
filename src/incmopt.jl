
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






