


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


# computes new candidate vectors for each chain
# accepts/rejects that vector on each chain, according to some rule
# *) computes N new parameter vectors
# *) applies a criterion to accept/reject any new params
# *) stores the result in chains
function computeNextIteration( algo::MAlgoRandom  )
  # here is the meat of your algorithm:
  # how to go from p(t) to p(t+1) ?
  for ch in 1:algo["N"]
  	# this updating rule can differ by chain!
  	algo.current_param[ch] = algo.current_param[ch] + randn()*shock_var
  end

    # check if we reached end of chain
	if algo.i == length(getEvals(algo.chains,1))
	    println("reached end of chain")
	    return true
	else
		# else update iteration on all chains
		updateIter!(algo.chains)

		if algo["mode"] == "serial"

			# map over 1:n evaluatceChainID
			v = map( x -> evaluateChainID(algo,x), 1:algo["N"])

		else

			# pmap over 1:n evaluateChainID
			v = pmap( x -> evaluateChainID(algo,x), 1:algo["N"])

		end
	end
	acceptReject!(algo)
	# append value to storage
	appendEval!(algo.chains,v)
end







