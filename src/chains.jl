# Types and methods related to the storing
# and updating of chains

type MChain

  n :: Int # number of chains
  i :: Int # current position in the chain
  L :: Int # total length of the chain

  evals   :: DataFrame   # DataFrame(id,iter,value)
  params  :: DataFrame
  moments :: DataFrame   # all previous params, indexed by iter

  function MChain(n,L,param_names,moment_names)
    new(n,1,L,DataFrame(),DataFrame(),DataFrame())
  end
end


function appendEval!(mc::MChain, value, moments, params)
  mc.i = mc.i +1
  mc.evals [mc.i]  = value
  mc.params[mc.i]  = collect(values(params))
  mc.moments[mc.i] = collect(values(moments))
end

