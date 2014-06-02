

type MCMChain

	id   :: Int
	iter :: Int
	p    :: Dict 	# current param value
	data :: DataFrame 	# DataFrame(id,iter,value)
	phist :: Dict 	# all previous params, indexed by iter

	function MCMChain(p,id)
		new(id,0,p,DataFrame(),[0 => p])
	end
end


type Moptim

	# Moptim setup
	initial_value    :: Dict 	# initial parameter value as a dict
	params_to_sample :: DataFrame 	
	objfunc          :: Function # objective function
	moments          :: DataFrame	# columns: data, sd
	moments_to_use   :: Array{ASCIIString,1} 	#names of moments to use
	shock_var        :: Float64 # variance of shock
	np_shock         :: Float64
	save_freq        :: Int
	N                :: Int 	# number of chains
	n_untempered     :: Int 	# number of chains untemprered
	maxiter          :: Int 	# maximum number of iterations
	run              :: Int 	# number of run
	i                :: Int 	# ?

	# cluster setup
	mode             :: ASCIIString	# {'serial','mpi'}

	# paths for I/O
	paths :: Dict

	prepared      :: Bool
	current_param :: Dict 	# current parameter value
	chains        :: Dict 	# collection of MCMChain's 

	# constructor
	function Moptim(initial_value,params_to_sample,objfunc,moments; moments_to_use=moments[:name],shock_var=1.0,np_shock=1.0,save_freq=25,N=3,n_untempered=1,mode="serial",paths=["chain"=> ".","lastparam"=>".","errorparam"=>".","wd"=>".","include_on_workers"=>"workers.jl","logdir"=>"./log"])

		iter = 0
		run = 0
		i   = 0

		# test format of moments_
		@assert haskey(moments.colindex,:name)
		@assert haskey(moments.colindex,:data)
		@assert haskey(moments.colindex,:sd)

		# names in params_to_sample_ must be members of key(initial_value_)
		@assert length(setdiff(collect(keys(initial_value)),params_to_sample[:name])) == 0

		# assign initial param to current param
		current_param = initial_value

		# setup chains
		chains = {i => MCMChain(current_param,i) for i = 1:N}

		if mode == "serial"
			prepared = true
		else
			prepared = false
		end


		return new(initial_value,params_to_sample,objfunc,moments,moments_to_use,shock_var,np_shock,save_freq,N,n_untempered,iter,run,i,mode,paths,prepared,current_param,chains)

	end	# constructor
end	#type

function show(io::IO,m::Moptim)
	print(io,"Moptim Object:\n")
	print(io,"==============\n\n")
	print(io,"Parameters to sample:\n")
	print(io,m.params_to_sample)
	print(io,"\nMoment Table:\n")
	print(io,m.moments)
	print(io,"\nMoment to use:\n")
	print(io,m.moments_to_use)
	print(io,"\nMode: $(m.mode)\n")
	print(io,"\nobjective function: $(m.objfunc)\n")
	print(io,"\nobjective function call: $(Expr(:call,m.objfunc,m.current_param,m.moments,m.moments_to_use))\n")
	if !m.prepared
		print(io,"\ncall MoptPrepare(m) to setup cluster\n")
	else 
		print(io,"\nNumber of chains: $(m.N)\n")
	end
	print(io,"END SHOW\n")
	print(io,"===========================\n")
end


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


# evaluate the objective function
# on all chains at their respective
# parameter values
function evaluateObjective(m::Moptim)

	if m.mode=="serial"
		# vals = eval(Expr(:call,m.objfunc,m.current_param,m.moments,m.moments_to_use))
		vals = map( x -> evalChain!(x,m) , values(m.chains))
	else
		if !m.prepared
			error("you must call MoptPrepare before evaluating the objective")
		end

		# pmap?
		vals = pmap( x -> evalChain!(x,m) , values(m.chains))
		# m.chains is a dict of different m.current_param. 
		# each different m.current_param represents a chain
	end
	return vals
end




# evalute chain
# calls the objective function with the 
# parameter vector currently stored in that chain
function evalChain!(c::MCMChain,m::Moptim)

	# call objective function of m
	x = eval(Expr(:call,m.objfunc,c.p,m.moments,m.moments_to_use))

	# store results in c
	c.data = rbind(c.data,DataFrame(id=c.id,iter=c.iter,value=x))
	c.phist[c.iter] = c.p

	return nothing
end

# update param on chain
function updateChain!(c::MCMChain,p::Dict)
	c.p = p
	c.iter += 1
	return nothing
end

# update param
function updateParam(m::Moptim)

	# for all chains in m.chains
	# compute a new guess

	# this is the actual MCMC algorithm
	# for each chain takes current and last value
	# and accepts/rejects it

end


# MCMChain methods
function summaryChain(c::MCMChain)

	# compute summary stats of chain

end





# define Test objective function
# function Testobj(x::Dict,mom::DataFrame,whichmom::Array{ASCIIString,1})

# 	mm = copy(mom)
# 	nm = names(mm)

# 	# add the model moments
# 	mm = cbind(mm,[x["a"] + x["b"] + rand() for i=1:nrow(mm)])
# 	names!(mm,[nm, :model])

# 	# subset to required moments only
# 	mm = mm[findin(mm[:name],whichmom),:]

# 	# compute distance
# 	v = sum((mm[:data]-mm[:model])^2)
# end
