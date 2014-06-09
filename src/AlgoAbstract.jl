# the algorithm should be the object
# created by the user. This is where the 
# argument should be set.

abstract type MAlgo
#   shock_var        :: Float64 # variance of shock
#   np_shock         :: Float64
#   save_freq        :: Int
#   N                :: Int   # number of chains
#   n_untempered     :: Int   # number of chains untemprered
#   maxiter          :: Int   # maximum number of iterations
#   run              :: Int   # number of run
#   i                :: Int   # ?

#   # cluster setup
#   mode             :: ASCIIString # {'serial','mpi'}

#   # paths for I/O
#   paths :: Dict

#   prepared      :: Bool
#   current_param :: Dict   # current parameter value
#   chains        :: MCMChain   # object of type MCMChain
    # paths=[
    #   "chain"=> ".",
    #   "lastparam"=>".",
    #   "errorparam"=>".",
    #   "wd"=>".",
    #   "include_on_workers"=>"workers.jl",
    #   "logdir"=>"./log"])
# end

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

type MAlgoRandom <: MAlgo
  opts :: Dict
  i :: Int

  function ChainAlgoRandom()
    return new(Dict(),0)
  end
end

function computeNewGuess( algo::ChainAlgoRandom  )
  algo.i = algo.i +1
end


