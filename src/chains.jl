
# the goal here is to have a convenient way of 
# storing the different evaluations.
# That includes the realizations of the objective
# function, the moments, and the value of the parameters
# that went in the evaluation.
# It might also include additional information
#
# so usually we'll have to refer to a chain number, an iteration, and an object
# which might be either a moment, a parameter, or an output.
#
# Because the chains might grow over time and because storing a in a DataFrame 
# is not always easy, we are going to use vectors for each.


# defining an abstract chain type in case 
# some algorithm need additional information
# on top of (eval, moment, parameter)
abstract AbstractChain


# the default chain type
# we create a dictionary with arrays
# for each parameters
type Chain
  i::Int             # current index
  evals     ::DataArray   # DataArray of evaluations (can hold NA)
  accept    ::DataArray   # DataArray of accept/reject(can hold NA)
  parameters::Dict   # dictionary of arrays(L,1), 1 for each parameter
  moments   ::Dict      # dictionary of DataArrays(L,1), 1 for each moment

  function Chain(MProb,L)
    evals      = @data([0.0 for i = 1:L])
    accept     = @data([false for i = 1:L])
    parameters = { x => zeros(L) for x in ps_names(MProb) }
    moments = { x => @data([0.0 for i = 1:L]) for x in ms_names(MProb) }
    return new(1,evals,accept,parameters,moments)
  end
end

# append an evaluation, the parameters and the moment to the stored data
# TODO let objfunc return a dict(eval,moments,params)
function appendEval!(chain::Chain, vals::Dict)
  chain.evals[chain.i] = vals["value"]
  for (k,v) in vals["moments"]
    chain.moments[k][chain.i] = v
  end
  for (k,v) in vals["params"]
    chain.params[k][chain.i] = v
  end
  return nothing
end


## MULTIPLE CHAINS
## ===============

# Stores multilpe chains
type MChain
  n :: Int # number of chains
  chains :: Array

  function MChain(n,MProb)
    chains = [  Chain(Mprob) for i in 1:n ]
    return new(n,chains)
  end
end


# evaluating the objective
function updateChain!(chain::Chain,m::MProb,p::Dict)

    # update counter on chain
    chain.i += 1

    # evaluate objective function
    v = eval(Expr(:call,m.objfunc,p,m.moments,m.moments_subset))

    # append to chain
    appendEval!(chain,v)

end




## -------------- OLD ----------------------------


type MCMChain

  evals   :: DataFrame   # DataFrame(id,iter,value)
  params  :: DataFrame 	 # DataFrame(id,iter,p1,p2,...)
  moments :: DataFrame 	 # DataFrame(id,iter,m1,m2,...)  

  function MCMChain(initial_value,N,moments)
  	devals   = DataFrame(id = [1:N],iter = [0 for i = 1:N],value = [0.0 for i = 1:N])
  	dparams  = DataFrame(id = [1:N],iter = [0 for i = 1:N])
  	dmoments = DataFrame(id = [1:N],iter = [0 for i = 1:N])

  	nms  = collect(keys(initial_value))
  	vals = collect(values(initial_value))

  	# build params dataFrame: one column for each
  	# elt of initial_value
  	for j in 1:length(initial_value)
  		dparams = cbind(dparams,[vals[j] for i=1:N])
  	end
  	names!(dparams, [:id,:iter,[symbol(nms[i]) for i=1:length(nms)]])

  	for j in 1:nrow(moments)
  		dmoments = cbind(dmoments,[0.0 for i=1:N])
  	end
  	names!(dmoments, [:id,:iter,Symbol[moments[i,:name] for i=1:nrow(moments)]])

    new(devals,dparams,dmoments)
  end
end


 # show methods
function showEvals(c::MCMChain)
	println(c.evals) 
	return nothing
end
function showParams(c::MCMChain)
	println(c.params) 
	return nothing
end
function showMoments(c::MCMChain)
	println(c.moments) 
	return nothing
end

# append methods
# take output of evaluateObjective and append to data in m.chains
function append!(c::MCMChain,v::Dict)

	N = length(v)

	# 1. append values
	c.evals = rbind(c.evals, c.evals[(end-N):end,:])
	c.evals[(end-N):end,:iter] .+= 1
	curriter = c.evals[end,iter]

	c.evals[(end-length(v)):end,:value] = [v[i]["values"] for i=1:N]

	# 2. append params
	tmp = c.params[(end-N):end,:]	# copy the last iteration
	tmp[:,:iter] .+= 1		# update iteration
	# get param names
	pnames = collect(keys(v[1]["param"]))
	for id =1:N
		for p in pnames
			tmp[tmp[:id]==id,symbol(p)] = v[id]["params"][pnames]
		end
	end
	c.params = rbind(c.params,tmp)	# append


	# 3. append moments
	tmp = c.moments[(end-N):end,:]
	tmp[:,:iter] .+= 1
	# get moments names
	mnames = v[1]["moments"][:name]
	for id =1:N
		for p in mnames
			tmp[tmp[:id]==id,symbol(p)] = v[id]["moments"][v[id]["moments"][:,:name]==p , :model]
		end
	end
	c.moments = rbind(c.moments,tmp)
	return nothing

end



