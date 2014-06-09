# the algorithm should be the object
# created by the user. This is where the 
# argument should be set.


# here is a list of potential options we want to set on each algo:

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
#   path :: ASCIIString # "where/to/save/this"

#   prepared      :: Bool
#   current_param :: Dict   # current parameter value
#   chains        :: MCMChain   # object of type MCMChain
# end




# define an abstract type and set/get for it
abstract MAlgo

# getter and setters for Algo
function getindex(algo::MAlgo, key)
  return(algo.opts[key])
end

function setindex!(algo::MAlgo, val,key)
  algo.opts[key] = val
end

function computeNewGuess( algo::MAlgo  )
  error("computeNewGuess not implemented for this algorithm")
end


# An implementation example
# -------------------------

type MAlgoRandom <: MAlgo
  m    :: MProb # an MProb
  opts :: Dict	# list of options
  i    :: Int 	# iteration
  current_param :: Dict
  chains :: MChain

  function MAlgoRandom(m::MProb,opts=["N"=>3,"shock_var"=>1.0,mode="serial","maxiter"=100],current_param=["a"=>1.1,"b" => 1.3])

  	# create chains
  	chains = MChain(opts["N"],m,opts["maxiter"])

    return new(m,opts,0,current_param,chains)
  end
end

function computeNewGuess( algo::MAlgoRandom  )
  algo.i = algo.i +1	# update iteration
  algo.current_param = algo.current_param + shock_var
end



# TODO: make a method of MAlgo
# evaluating the objective
# and appendEval
function updateChain!(chain::Chain,m::MProb,p::Dict)

    if chain.i == length(chain.evals)
        println("reached end of chain")
        return true
    else
        # update counter on chain
        chain.i += 1

        # evaluate objective function
        v = eval(Expr(:call,m.objfunc,p,m.moments,m.moments_subset))

        # append to chain
        appendEval!(chain,v)

    end


end





