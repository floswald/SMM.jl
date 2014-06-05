# the algorithm should be the object
# created by the user. This is where the 
# argument should be set.

abstract ChainAlgoAbstract

# getter and setters for Algo
function getindex(algo::ChainAlgoAbstract, key)
  return(algo.opts[key])
end

function setindex!(algo::ChainAlgoAbstract, val,key)
  algo.opts[key] = val
end

function computeNewGuess( algo::ChainAlgoAbstract  )
  error("computeNewGuess not implemented for this algorithm")
end


# An implementation example
# -------------------------

type ChainAlgoRandom <: ChainAlgoAbstract
  opts :: Dict
  i :: Int

  function ChainAlgoRandom()
    return new(Dict(),0)
  end
end

function computeNewGuess( algo::ChainAlgoRandom  )
  algo.i = algo.i +1
end


