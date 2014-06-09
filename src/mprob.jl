# the type MProb describes the optimization
# problem, it knows about the parameters
# to estimate, their names and bounds, 
# but also about the moments and their names.
type MProb

  # Moptim setup
  initial_value    :: Dict  # initial parameter value as a dict
  params           :: Dict  # a dictionary of upper and lower bound   
  moments          :: Dict  # a dictionary of moments to track
  param_init       :: Dict
  objfunc          :: Function # objective function

  # constructor
  function Moptim(
    initial_value,
    params_to_sample,
    objfunc,moments; 
    moments_to_use=moments[:name],
    shock_var=1.0,
    np_shock=1.0,
    save_freq=25,
    N=3,
    n_untempered=1,
    mode="serial",
    paths=[
      "chain"=> ".",
      "lastparam"=>".",
      "errorparam"=>".",
      "wd"=>".",
      "include_on_workers"=>"workers.jl",
      "logdir"=>"./log"])

    iter = 0
    run = 0
    i   = 0

    # test format of moments_
    @assert haskey(moments.colindex,:name)
    @assert haskey(moments.colindex,:data)
    @assert haskey(moments.colindex,:sd)

    # names in params_to_sample_ must be members of key(initial_value_)
    @assert length(setdiff(collect(keys(initial_value)),params_to_sample[:name])) == 0

    # assign initial param to current param for each chain
    current_param = {i => initial_value for i=1:N}

    # setup chains
    chains = MCMChain(initial_value,N,moments)

    if mode == "serial"
      prepared = true
    else
      prepared = false
    end


    return new(initial_value,params_to_sample,objfunc,moments,moments_to_use,shock_var,np_shock,save_freq,N,n_untempered,iter,run,i,mode,paths,prepared,current_param,chains)

  end # constructor
end #type

# returns the list of paramaters to sample
function ps_names(mprob::MProb)
  return(keys(mprob.params))
end

function ms_names(mprob::MProb)
  return(keys(mprob.moments))
end




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