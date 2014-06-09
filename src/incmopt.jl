
# define a Test objective function
function Testobj(x::Dict,mom::Dict,whichmom::Array{ASCIIString,1})

	t0 = time()
	mm = DataFrame(name=collect(keys(mom)),data=map(x -> x[1],values(mom)),model=[x["a"] + x["b"] + i for i=1:length(mom)])

	# get model moments as a dict
	mdict = {i => mm[findin(mm[:name],[i]),:model][1] for i in collect(keys(mom)) }

	# subset to required moments only
	mm = mm[findin(mm[:name],whichmom),:]

	# compute distance
	v = sum((mm[:data] - mm[:model]).^2)

	# status
	status = 1

	# time out
	t0 = time() - t0

	# return a dict
	ret = ["value" => v, "params" => x, "time" => t0, "status" => status, "moments" => mdict]
	return ret

end



# # call this function to set up
# # the cluster. 
# # maps choice of "mode" into an action
# function MoptPrepare!(m::Moptim)

# 	if m.prepared
# 		println("you're good to go")
# 		return nothing
# 	else

# 		# setup cluster
# 		@everywhere include(m.include_on_workers)

# 		# set m.N = # of workers

# 		# send data to workers

# 		# set m.prepared = true
# 	end

# end


# function evaluateObjective(m::Moptim)

# 	if m.mode =="serial"
# 		v = map(x -> evaluateChainID(x,m), 1:m.N )
# 	else
# 		v = pmap(x -> evaluateChainID(x,m), 1:m.N )

# 	end
# 	return v

# end

# function evaluateChainID(m::Moptim,i::Int)

# 	# eval chain i with param[i]
# 	x = eval(Expr(:call,m.objfunc,m.current_param[i],m.moments,m.moments_to_use))
# 	return x

# end

# function collectChains(m::Moptim)



# end






