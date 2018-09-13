
"""
# Minimisation Problem: `MProb`

A moment **minimsation** problem is defined by an objective function that
depends on a vector of unknown parameters `params_to_sample`, and a set
of datamoments `moments`. 

## Fields:

* `initial_value`: initial parameter value as a dict
* `params_to_sample`: OrderedDict with lower and upper bounds
* `objfunc`: objective function
* `objfunc_opts`: options passed to the objective function, e.g. printlevel
* `moments`: a dictionary of data moments to track

## Example:

    ```julia
    pb    = Dict(p1" => [0.2,-2,2] , "p2" => [-0.2,-2,2] )
    moms  = DataFrame(name=["mu2","mu1"],value=[0.0,0.0],weight=rand(2))
    m     = MProb() 
    addSampledParam!(m,pb) 
    addMoment!(m,moms) 
    addEvalFunc!(m.MOpt.objfunc_norm)
    ```
    
"""
mutable struct MProb

  # setup
  initial_value       :: OrderedDict           # initial parameter value as a dict
  params_to_sample    :: OrderedDict           # Dict with lower and upper bounds
  objfunc             :: Function       # objective function
  objfunc_opts        :: Dict           # options passed to the objective function, e.g. printlevel
  moments             :: OrderedDict           # a dictionary of data moments to track

  # very simple constructor
  function MProb()
    this = new()
    this.initial_value       = OrderedDict()
    this.params_to_sample    = OrderedDict()
    function def(x)
      @info "default objfunc, returns input"
      x
    end
    this.objfunc             = def
    this.objfunc_opts        = Dict()
    this.moments             = OrderedDict()
    return(this)
  end

end #type


# -------------------- ADDING PARAMS --------------------
"""
Add initial parameter values to an `MProb` minimsation problem.
"""
function addParam!(m::MProb,name::String,init)
  m.initial_value[Symbol(name)] = init
end

"""
Add initial parameter values to an `MProb` minimsation problem.

### Arguments:

* `p`: A Dict with (String,Number) pairs
"""
function addParam!(m::MProb,p::Union{Dict,OrderedDict})
  for k in keys(p)
    m.initial_value[Symbol(k)] = p[k]
  end
end


"""
Add parameters to be sampled to an `MProb`.
"""
function addSampledParam!(m::MProb,name::Any,init::Any, lb::Any, ub::Any)
  @assert ub>lb
  m.initial_value[Symbol(name)] = init 
  m.params_to_sample[Symbol(name)] = Dict( :lb => lb , :ub => ub)
  return m
end

"""
Add parameters to be sampled to an `MProb`.

`d`: a Dict with a triple (init,lb,ub) as value for each key.
"""
function addSampledParam!(m::MProb,d=OrderedDict{Any,Array{Any,1}})
  for k in keys(d)
    addSampledParam!(m,Symbol(k),d[k][1],d[k][2],d[k][3])
  end
  return m
end

# -------------------- ADDING MOMENTS --------------------

"""
Add Moments to an `MProb`.

add a single moment to the mprob.

`name`: the name of the moment as a Symbol
`value`: value of the moment
`weight`: weight in the objective function
"""
function addMoment!(m::MProb,name::Symbol,value,weight)
  m.moments[Symbol(name)] = Dict( :value => value , :weight => weight )
  return m 
end


addMoment!(m::MProb,name::String,value,weight) = addMoment!(m,Symbol(name),value,weight)
addMoment!(m::MProb,name::String,value) = addMoment!(m,Symbol(name),value,1.0)
addMoment!(m::MProb,name::Symbol,value) = addMoment!(m,(name),value,1.0)

function addMoment!(m::MProb,d::Dict)
  for k in keys(d)
    addMoment!(m,Symbol(k),d[k][:value],d[k][:weight])
  end
  return m 
end

function addMoment!(m::MProb,d::DataFrame)
  for i in 1:nrow(d)
    m = addMoment!(m,Symbol(d[i,:name]),d[i,:value],d[i,:weight])
  end
  return m 
end


# -------------------- ADDING OBJ FUNCTION --------------------
function addEvalFunc!(m::MProb,f::Function)
  m.objfunc = f
end

# evalute objective function
function evaluateObjective(m::MProb,p::Union{Dict,OrderedDict})
    ev = Eval(m,p)
    try
       # ev = eval(Expr(:call,m.objfunc,ev))
      ev = m.objfunc(ev)
    catch ex
      @warn "caught exception $ex"
      ev.status = -2
    end
    gc()
    return ev
end

function evaluateObjective(m::MProb,ev)
    # catching errors
    try
       # ev = eval(Expr(:call,m.objfunc,ev))
       ev = m.objfunc(ev)
    catch ex
      @warn "caught exception $ex"
      ev.status = -2
    end
    gc()
    return ev
end

# -------------------- GETTERS --------------------

"""
Get all parameter names from the `MProb`
"""
function ps_names(mprob::MProb)
  return(keys(mprob.initial_value))
end

"""
Get sampled parameter names from the `MProb`
"""
function ps2s_names(mprob::MProb)
  return convert(Array{Symbol,1},[k for k in keys(mprob.params_to_sample)])
end

"""
Get the name of moments
"""
function ms_names(mprob::MProb)
  return(keys(mprob.moments))
end

function show(io::IO,m::MProb)

  print(io,"MProb Object:\n")
  print(io,"==============\n\n")
  print(io,"Parameters to sample:\n")
  print(io,m.params_to_sample)
  print(io,"\nMoment Table:\n")
  print(io,m.moments)
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
