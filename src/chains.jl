# Types and methods related to the storing
# and updating of chains

type MCMChain

  n :: Int # number of chains

  evals   :: DataFrame   # DataFrame(id,iter,value)
  params  :: DataFrame
  moments :: DataFrame   # all previous params, indexed by iter

  function MCMChain(p,id)
    new(id,0,p,DataFrame(),[0 => p])
  end
end