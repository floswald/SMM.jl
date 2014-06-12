# the type MProb describes the optimization
# problem, it knows about the parameters
# to estimate, their names and bounds, 
# but also about the moments and their names.
type MProb

  # setup
  initial_value    :: Dict  # initial parameter value as a dict
  # params_to_sample :: Dict{ASCIIString,Array{Float64,1}}  # a dictionary of upper and lower bound for params we estimate (others are fixed)
  params_to_sample :: Dict
  objfunc          :: Function # objective function
  moments          :: Dict  # a dictionary of moments to track
  moments_subset   :: Array{ASCIIString}  # an array of moment names to subset objective funciton

  # constructor
  function MProb(
    initial_value,
    params_to_sample,
    objfunc,moments; 
    moments_subset=collect(keys(moments)))

    # assert that moments_subset is a subset of moment keys
    @assert issubset(moments_subset, keys(moments))

    # assert that moments has two entries for each moment: value and sd
    @assert all(map(x -> length(x),values(moments)) .== 2)

    # assert that params_to_sample are subset of initial_value
    @assert issubset(keys(initial_value),keys(params_to_sample))

    # assert that params_to_sample has valid bounds: a vector with 2 increasing entries
    @assert all(map(x -> x[2]-x[1],values(params_to_sample)) .> 0)

    return new(initial_value,params_to_sample,objfunc,moments,moments_subset)

  end # constructor
end #type

# returns the list of paramaters to sample
function ps_names(mprob::MProb)
  return(keys(mprob.initial_value))
end

function ms_names(mprob::MProb)
  return(keys(mprob.moments))
end




function show(io::IO,m::MProb)

  mdf = DataFrame(moment=collect(keys(m.moments)),value=map(x -> x[1], values(m.moments)),sd=map(x -> x[2], values(m.moments)))
  pdf = DataFrame(param=collect(keys(m.params_to_sample)),lb=map(x -> x[1], values(m.params_to_sample)),ub=map(x -> x[2], values(m.params_to_sample)))

  print(io,"MProb Object:\n")
  print(io,"==============\n\n")
  print(io,"Parameters to sample:\n")
  print(io,pdf)
  print(io,"\nMoment Table:\n")
  print(io,mdf)
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