# the algorithm should be the object
# created by the user. This is where the 
# argument should be set.


# here is a list of potential options we may want to set on each algo:

#   shock_var        :: Float64 # variance of shock
#   np_shock         :: Float64
#   save_freq        :: Int
#   N                :: Int   # number of chains
#   n_untempered     :: Int   # number of chains untemprered
#   maxiter          :: Int   # maximum number of iterations
#   run              :: Int   # number of run

#   # cluster setup
#   mode             :: ASCIIString # {'serial','mpi'}

#   # paths for I/O
#   path :: ASCIIString # "where/to/save/this"

#   prepared      :: Bool
#   current_param :: Dict   # current parameter value
#   chains        :: MChain   # object of type MChain
# end




# define an abstract type and set/get for it
abstract MAlgo

# getter and setters for Algo.opts
function getindex(algo::MAlgo, key)
  return(algo.opts[key])
end

function setindex!(algo::MAlgo, val,key)
  algo.opts[key] = val
end

function computeNewGuesses( algo::MAlgo  )
  error("computeNewGuesses not implemented for this algorithm")
end


# An implementation example
# -------------------------

# user defines their MCMC algorithm "Random"
# requires:  
# 1) type declaration 
# 2) method computeNewGuess

type MAlgoRandom <: MAlgo
  m             :: MProb # an MProb
  opts          :: Dict	# list of options
  i             :: Int 	# iteration
  current_param :: Array{Dict,1}  # current param value: one Dict for each chain
  chains        :: MChain 	# collection of Chains: if N==1, length(chains) = 1

  function MAlgoRandom(m::MProb,opts=["N"=>3,"shock_var"=>1.0,"mode"=>"serial","maxiter"=>100])

  	# create chains
  	chains = MChain(opts["N"],m,opts["maxiter"])
  	cpar = [ m.initial_value for i=1:opts["N"] ] 

    return new(m,opts,0,cpar,chains)
  end
end


# computes new guesses for each chain
# i.e. computes N new parameter vectors
function computeNewGuesses( algo::MAlgoRandom  )
  algo.i += 1	# update iteration on algorithm
  # here is the meat of your algorithm:
  # how to go from p(t) to p(t+1) ?
  for ch in 1:algo["N"]
  	# this updating rule can differ by chain!
  	algo.current_param[ch] = algo.current_param[ch] + shock_var
  end
end



# updates all chains in algo
# wrapper to evaluateChainID
# calls map(evaluateChainID) or pmap(evaluateChainID) depending on opts
function updateChains!(algo::MAlgo)

	# check if we reached end of chain
    if algo.chains[1].i == length(chain.evals)
        println("reached end of chain")
        return true
    else
    	if algo["mode"] == "serial"

    		# map over 1:n evaluatceChainID
    		v = map( x -> evaluateChainID(algo,x), 1:algo["N"])

    	else

    		# pmap over 1:n evaluateChainID
    		v = pmap( x -> evaluateChainID(algo,x), 1:algo["N"])

    	end
    end
    # append value to storage
    appendEval!(algo.chains,v)
end
    

# evalute chain number i
# with param vector number i
function evaluateChainID(algo::MAlgo,i::Int)

	# update iteration on chain i
	algo.chains[i] =+ 1

	# eval chain i with param[i]
	x = eval(Expr(:call,algo.m.objfunc,algo.current_param[i],algo.m.moments,algo.m.moments_subset))
	return x

end




