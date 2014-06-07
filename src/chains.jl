# Types and methods related to the storing
# and updating of chains
type MChain

  n :: Int # number of chains
  i :: Int # current position in the chain
  L :: Int # total length of the chain

    function MChain(n,L,param_names,moment_names)
    new(n,1,L,DataFrame(),DataFrame(),DataFrame())
  end
end

function appendEval!(mc::MChain, value, moments, params)
  mc.i = mc.i +1
  mc.evals [mc.i]  = value
  mc.params[mc.i]  = collect(values(params))	# that does not append row-wise
  mc.moments[mc.i] = collect(values(moments))
end


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



