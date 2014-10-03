
export addMoment


#'.. py:class:: MProb
#'
#'   the type MProb describes the optimization
#'   problem, it knows about the parameters
#'   to estimate, their names and bounds, 
#'   but also about the moments and their names.
type MProb

  # setup
  initial_value       :: Dict           # initial parameter value as a dict
  params_to_sample    :: Dict           # Dict with lower and upper bounds
  objfunc             :: Function       # objective function
  objfunc_opts        :: Dict           # options passed to the objective function, e.g. printlevel
  moments             :: Dict           # a dictionary of moments to track
  moments_subset      :: Array{Symbol}  # an array of moment names to subset objective funciton

  # very simple constructor
  function MProb()
    this = new()
    this.initial_value       = Dict()
    this.params_to_sample       = Dict()
    this.objfunc             = x -> x
    this.objfunc_opts        = Dict()
    this.moments             = Dict()
    this.moments_subset      = []
    return(this)
  end

end #type


# -------------------- ADDING PARAMS --------------------
function addParam!(m::MProb,name::ASCIIString,init)
  m.initial_value[symbol(name)] = init
  return m 
end

function addSampledParam!(m::MProb,name::Any, init::Any, lb::Any, ub::Any)
  @assert ub>lb
  m.initial_value[symbol(name)] = init 
  m.params_to_sample[symbol(name)] = [ :lb => lb , :ub => ub]
  return m
end


# -------------------- ADDING MOMENTS --------------------

function addMoment(m::MProb,name::Symbol,value,weight)
  m.moments[symbol(name)] = [ :value => value , :weight => weight ]
  return m 
end
addMoment(m::MProb,name::ASCIIString,value,weight) = addMoment(m,symbol(name),value,weight)
addMoment(m::MProb,name::ASCIIString,value) = addMoment(m,symbol(name),value,1.0)
addMoment(m::MProb,name::Symbol,value) = addMoment(m,(name),value,1.0)

function addMoment(m::MProb,d::Dict)
  for k in keys(d)
    addMoment(m,symbol(k),d[k][:value],d[k][:weight])
  end
  return m 
end

function addMoment(m::MProb,d::DataFrame)
  for (i in 1:nrow(d))
    m = addMoment(m,symbol(d[i,:name]),d[i,:value],d[i,:weight])
  end
  return m 
end

# -------------------- GETTERS --------------------

#'.. py:function:: ps_names
#'
#'   returns the list of paramaters to sample
function ps_names(mprob::MProb)
  return(keys(mprob.initial_value))
end
function ps2s_names(mprob::MProb)
  return(keys(mprob.params_to_sample))
end

function ms_names(mprob::MProb)
  return(keys(mprob.moments))
end

function show(io::IO,m::MProb)

  print(io,"MProb Object:\n")
  print(io,"==============\n\n")
  print(io,"Parameters to sample:\n")
  print(io,m.params_to_sample_df)
  print(io,"\nMoment Table:\n")
  print(io,m.moments)
  print(io,"Moment to use:\n")
  print(io,m.moments_subset)
  print(io,"\n")
  print(io,"\nobjective function: $(m.objfunc)\n\n")
  # print(io,"\nobjective function call: $(Expr(:call,m.objfunc,m.current_param,m.moments,m.moments_to_use))\n")
  # if !m.prepared
  #   print(io,"\ncall MoptPrepare(m) to setup cluster\n")
  # else 
  #   print(io,"\nNumber of chains: $(m.N)\n")
  # end
  print(io,"\nEND SHOW MProb\n")
  print(io,"===========================\n")
end
