
# The BGP MCMC Algorithm: Likelihood-Free Parallel Tempering
# ==========================================================
#
# http://link.springer.com/article/10.1007%2Fs11222-012-9328-6
# http://fr.arxiv.org/abs/1108.3423
#
# Baragatti, Grimaud and Pommeret (BGP)
# 
# Approximate Bayesian Computational (ABC) methods (or likelihood-free methods) have appeared in the past fifteen years as useful methods to perform Bayesian analyses when the likelihood is analytically or computationally intractable. Several ABC methods have been proposed: Monte Carlo Markov Chains (MCMC) methods have been developped by Marjoramet al. (2003) and by Bortotet al. (2007) for instance, and sequential methods have been proposed among others by Sissonet al. (2007), Beaumont et al. (2009) and Del Moral et al. (2009). Until now, while ABC-MCMC methods remain the reference, sequential ABC methods have appeared to outperforms them (see for example McKinley et al. (2009) or Sisson et al. (2007)). In this paper a new algorithm combining population-based MCMC methods with ABC requirements is proposed, using an analogy with the Parallel Tempering algorithm (Geyer, 1991). Performances are compared with existing ABC algorithms on simulations and on a real example.


# Define a Chain Type for BGP
type BGPChain <: AbstractChain
  id::Int             # chain id
  i::Int              # current index
  infos      ::DataFrame   # DataFrameionary of arrays(L,1) with eval, ACC and others
  parameters ::DataFrame   # DataFrameionary of arrays(L,1), 1 for each parameter
  moments    ::DataFrame   # DataFrameionary of DataArrays(L,1), 1 for each moment
  accept_tol ::Float64     # acceptance tolerance: new is not a "substantial" improvement over old, don't accept
  dist_tol   ::Float64 # percentage value distance from current chain i that is considered "close" enough. i.e. close = ((val_i - val_j)/val_j < dist_tol )
  jump_prob  ::Float64 # probability of swapping with "close" chain

  params_nms ::Array{Symbol,1}	# names of parameters (i.e. exclusive of "id" or "iter", etc)
  moments_nms::Array{Symbol,1}	# names of moments

  # TODO need either of those not both
  # the paper uses tempering to set up the kernel,
  # tibo uses shock_sd to amplify the shocks. expect small difference.
  tempering  ::Float64 # tempering in update probability
  shock_sd   ::Float64 # sd of shock to 

  function BGPChain(id,MProb,L,temp,shock,accept_tol,dist_tol,jump_prob)
    infos      = DataFrame(chain_id = [id for i=1:L], iter=1:L, evals = zeros(Float64,L), accept = zeros(Bool,L), status = zeros(Int,L), exchanged_with=zeros(Int,L),prob=zeros(Float64,L),ratio_old_new=zeros(Float64,L),accept_rate=zeros(Float64,L),shock_sd = [shock,zeros(Float64,L-1)],eval_time=zeros(Float64,L))
    parameters = cbind(DataFrame(chain_id = [id for i=1:L], iter=1:L), convert(DataFrame,zeros(L,length(ps_names(MProb)))))
    moments = cbind(DataFrame(chain_id = [id for i=1:L], iter=1:L), convert(DataFrame,zeros(L,length(ms_names(MProb)))))
    par_nms = sort(Symbol[ symbol(x) for x in ps_names(MProb) ])
    mom_nms = sort(Symbol[ symbol(x) for x in ms_names(MProb) ])
    names!(parameters,[:chain_id,:iter, par_nms])
    names!(moments   ,[:chain_id,:iter, mom_nms])
    # infos      = { "evals" => @data([0.0 for i = 1:L]) , "accept" => @data([false for i = 1:L]), "status" => [0 for i = 1:L], "exchanged_with" => [0 for i = 1:L]}
    # parameters = { x => zeros(L) for x in ps_names(MProb) }
    # moments    = { x => @data([0.0 for i = 1:L]) for x in ms_names(MProb) }
    return new(id,0,infos,parameters,moments,accept_tol,dist_tol,jump_prob,par_nms,mom_nms,temp,shock)
  end
end



type MAlgoBGP <: MAlgo
  m               :: MProb # an MProb
  opts            :: Dict	# list of options
  i               :: Int 	# iteration
  current_param   :: Array{Dict,1}  # current param value: one Dict for each chain
  candidate_param :: Array{Dict,1}  # dict of candidate parameters: if rejected, go back to current
  MChains         :: Array{BGPChain,1} 	# collection of Chains: if N==1, length(chains) = 1

  function MAlgoBGP(m::MProb,opts=["N"=>3,"min_shock_sd"=>0.1,"max_shock_sd"=>1.0,"mode"=>"serial","maxiter"=>100,"maxtemp"=> 100])

  	# temperatures for each chain
	temps = linspace(1.0,opts["maxtemp"],opts["N"])
  	# acceptance tolerance within each chain. condition: abs(old - new)/abs(old) > tol
	acctol = linspace(opts["min_accept_tol"],opts["max_accept_tol"],opts["N"])
  	# shock standard deviations for each chain
	shocksd = linspace(opts["min_shock_sd"],opts["max_shock_sd"],opts["N"])
  	# acceptance tolerance for cross chain jumps. condition: abs(val(1) - val(2)) < tol
	disttol = linspace(opts["min_disttol"],opts["max_disttol"],opts["N"])
  	# acceptance tolerance for cross chain jumps. condition: abs(val(1) - val(2)) < tol
	jump_prob = linspace(opts["min_jump_prob"],opts["max_jump_prob"],opts["N"])
  	# create chains
  	chains = [BGPChain(i,m,opts["maxiter"],temps[i],shocksd[i],acctol[i],disttol[i],jump_prob[i]) for i=1:opts["N"] ]
  	# current param values
  	cpar = [ deepcopy(m.initial_value) for i=1:opts["N"] ] 
  	# candidate param values
  	cpar0 = [ deepcopy(m.initial_value) for i=1:opts["N"] ] 
  	# jump register
  	

    return new(m,opts,0,cpar,cpar0,chains)
  end
end


function resetCurrentParam2initial(algo::MAlgoBGP)
	for j in 1:length(algo.current_param)
		algo.current_param[j] = algo.m.initial_value
	end
end



# computes new candidate vectors for each chain
# accepts/rejects that vector on each chain, according to some rule
# *) computes N new parameter vectors
# *) applies a criterion to accept/reject any new params
# *) stores the result in chains
function computeNextIteration!( algo::MAlgoBGP  )
    # here is the meat of your algorithm:
    # how to go from p(t) to p(t+1) ?


    # check if we reached end of chain
	if algo.i == algo["maxiter"]
	    println("reached end of chain. goodbye.")
	    return true
	else
		# else update iteration count on all chains
		updateIterChain!(algo.MChains)


	
		@assert algo.i == algo.MChains[1].i

		# New Candidates
		# --------------

		if algo.i > 1
			# MVN = getParamKernel(algo)	# returns a MvNormal object
			MVN = getParamCovariance(algo)	# returns a Cov matrix
			getNewCandidates!(algo,MVN)
		end

		# evaluate objective on all chains
		# --------------------------------
		if algo["mode"] == "serial"
			v = map( x -> evaluateObjective(algo,x), 1:algo["N"])
		else
			v = pmap( x -> evaluateObjective(algo,x), 1:algo["N"])
		end

		# notice that at this point, v[chain_id]["param"] is the candidate param vector
		# we append v to the chain
		# so if candidate is rejected, go back into v[chain_id]["param"] and reset to previous param


		# Part 1) LOCAL MOVES ABC-MCMC for i={1,...,N}. accept/reject
		# -----------------------------------------------------------

		localMovesMCMC!(algo,v)

		# # Part 2) EXCHANGE MOVES 
		# # ----------------------
		# if algo["N"] >1
		# 	exchangeMoves!(algo)
		# end

		# Part 3) update sampling variances
	end
end





# notice: higher tempering draws candiates further spread out,
# but accepts lower function values with lower proability
function localMovesMCMC!(algo::MAlgoBGP,v::Array{Dict{ASCIIString,Any},1})
	for ch in 1:algo["N"]
		if algo.i == 1
			# accept all
			ACC = true
	  		algo.current_param[ch] = deepcopy(algo.candidate_param[ch])
	  		status = 1
	  		prob = 1.0
			xnew = v[ch]["value"][1]
			xold = v[ch]["value"][1]
			algo.MChains[ch].infos[algo.i,:accept_rate] = 0.1
		    # append values to MChains at index ch
		    appendEval!(algo.MChains[ch],v[ch],ACC,status,prob)
		else
			xold = evals(algo.MChains[ch],algo.i-1)[1]
			xnew = v[ch]["value"][1]
			prob = minimum([1.0, exp(algo.MChains[ch].tempering *(xold - xnew))])
			prob = prob * (xnew < algo.MChains[ch].accept_tol)
			if isna(prob)
				prob = 0.0
				status = -1
				ACC = false
			elseif !isfinite(xold)
				prob = 1.0
				status = -2
				ACC = false
			else 
				if prob > rand()
					ACC = true
			  		algo.current_param[ch] = deepcopy(algo.candidate_param[ch] )
				else
					ACC = false
					v[ch]["value"] = xold # reset value in output of obj to previous value
					v[ch]["params"] = deepcopy(algo.current_param[ch])	# reset param in output of obj to previous value
					v[ch]["moments"] = df2dict(moments(algo.MChains[ch],algo.i-1))	# reset moments in output of obj to previous value
				end
				status = 1
			end 
			# append values to MChains at index ch
		    appendEval!(algo.MChains[ch],v[ch],ACC,status,prob)
		    # update sampling variances
		    algo.MChains[ch].infos[algo.i,:accept_rate] = 0.9 * algo.MChains[ch].infos[algo.i-1,:accept_rate] + 0.1 * ACC
		    algo.MChains[ch].shock_sd = algo.MChains[ch].shock_sd * (1+ 0.05*( 2*(algo.MChains[ch].infos[algo.i,:accept_rate]>0.234) -1) )
		    algo.MChains[ch].infos[algo.i,:shock_sd] = algo.MChains[ch].shock_sd
		    algo.MChains[ch].infos[algo.i,:ratio_old_new] = xold / xnew
		end
		if algo.i>1 && algo["N"] > 1 
			exchangeMoves!(algo,ch,xold)
		end
	end
	if mod(algo.i,100) == 0
		println(infos(algo.MChains,algo.i))
	end
end

function exchangeMoves!(algo::MAlgoBGP,ch::Int,oldval)

	# 1) find all chains with value +/- 10% of chain ch
	for ch2 in 1:algo["N"]
		dist = Float64[]
		idx = Int64[]
		if ch != ch2
			tmp = abs(evals(algo.MChains[ch2],algo.MChains[ch2].i)[1] - oldval) / abs(oldval)
			push!(dist,tmp)
			push!(idx,ch2)
		end
		close = idx[dist .< algo.MChains[ch].dist_tol]
		if length(close) >0
			# 2) with 5% probability exchange with a randomly chosen chain from close
			if rand() < algo.MChains[ch].jump_prob
				swapRows!(algo,(ch,sample(close)),algo.i)
			end
		end
	end

end




# 	if !algo["rings"]

# 		# Part 2a) EXCHANGE MOVES without rings
# 		# ======================
# 		# 1) N exchange moves are proposed
# 		# 2) N pairs are chosen uniformly among all possible pairs with replacment 
# 		# 3) exchange (z_i,theta_i) with (z_j,theta_j) if val_j < tol_i

# 		jump_pairs = sample(algo.Jump_register,algo["N"],replace=false)
# 		for pair in jump_pairs
# 			distance = abs(evals(algo.MChains[pair[1]],algo.MChains[pair[1]].i) - evals(algo.MChains[pair[2]],algo.MChains[pair[2]].i))
# 			# println("distance of pair $(pair) = $distance")
# 			if distance[1] < algo.MChains[pair[1]].jumptol
# 				swapRows!(algo,pair,algo.i)
# 			end
			
# 		end

# 	else
# 		# Part 2b) EXCHANGE MOVES with rings
# 		# ==================================
# 		# 1) N exchange moves are proposed
# 		# 2) a ring E_j is a partition of objective function space, i.e. 
# 		#    chain j \in E_1 iff V_j \in [0,E_1]
# 		#    chain j \in E_2 iff V_j \in [E_1,E_2]
# 		#    ...                             
# 		# 3) randomly choose N rings, with replacement. A ring with less than 2 chains is not admissible.                               
# 		# 5) from within that ring, form jump_pairs 
# 		# 6) exchange (z_i,theta_i) with (z_j,theta_j) if val_j < tol_i

# 		chain_in_rings = infos(algo.MChains,algo.i)[[:chain_id,:evals]]
# 		chain_in_rings = cbind(chain_in_rings,DataFrame(ring=findInterval(chain_in_rings[:evals],algo["rings"])))

# 		# drop rings with less than 2 chains
# 		hi = hist(chain_in_rings[:ring])
# 		keeps = findin(chain_in_rings[:ring],unique(sort(chain_in_rings[:ring]))[hi[2].>1])
# 		chain_in_rings = chain_in_rings[keeps,:]

# 		if nrow(chain_in_rings) > 0

# 			# draw pairs to exchange with 5% probability
# 			chain_in_rings = cbind(chain_in_rings,DataFrame(exchange=sample([true,false],WeightVec([0.05,0.95]),nrow(chain_in_rings))))
# 			chain_in_rings = chain_in_rings[chain_in_rings[:exchange] .== true,:]

# 		# 	# choose N rings with replacement
# 			rings = DataFrame(ringid=sample(chain_in_rings[:ring],algo["N"],replace=true), pair1 = [0 for i=1:algo["N"]], pair2 = [0 for i=1:algo["N"]])

# 			# for each entry of rings, randomly choose a pair of chains
# 			for ir = 1:nrow(rings)

# 					# @bp length(chain_in_rings[chain_in_rings[:ring] .== rings[ir,:ringid],:chain_id]) ==1


# 				pairs = sample(chain_in_rings[chain_in_rings[:ring] .== rings[ir,:ringid],:chain_id],2,replace=false)
# 				# add to rings
# 				rings[ir,:pair1] = pairs[1]
# 				rings[ir,:pair2] = pairs[2]

# 			end

# 			# for each ring, if val(j) < tol_i, swap i with j
# 			for irow in eachrow(rings)
# 				distance = abs(evals(algo.MChains[irow[:pair2]],algo.MChains[irow[:pair2]].i)[1]) < algo.MChains[irow[:pair1]].jumptol
# 				if distance
# 					swapRows!(algo,(irow[:pair1],irow[:pair2]),algo.i,irow[:ringid])
# 				end
# 			end
# 		end
# 	end
# end


function findInterval(x,vec::Array)

	out = zeros(Int,length(x))
	sort!(vec)
	i = 0

	for j in 1:length(x)
		out[j] = findfirst(vec .> x[j]) - 1
	end
	return out
end


function swapRows!(algo::MAlgoBGP,pair::(Int,Int),i::Int)

	# pars, moms and value from 1
	p1 = parameters(algo.MChains[pair[1]],i)
	m1 = moments(algo.MChains[pair[1]],i)
	v1 = evals(algo.MChains[pair[1]],i)

	# same for 2
	p2 = parameters(algo.MChains[pair[2]],i)
	m2 = moments(algo.MChains[pair[2]],i)
	v2 = evals(algo.MChains[pair[2]],i)

	# make a note in infos
	algo.MChains[pair[1]].infos[i,:exchanged_with] = pair[2]
	algo.MChains[pair[2]].infos[i,:exchanged_with] = pair[1]

	# swap
	algo.MChains[pair[1]].parameters[i,:] = p2
	algo.MChains[pair[2]].parameters[i,:] = p1
	algo.MChains[pair[1]].moments[i,:] = m2
	algo.MChains[pair[2]].moments[i,:] = m1
	algo.MChains[pair[1]].infos[i,:evals] = v2[1]
	algo.MChains[pair[2]].infos[i,:evals] = v1[1]

end


# get past parameter values to compute new candidates
# we draw new candidates for each chain from a joint normal
# that depends on parameter vectors on ALL chains
# pardf = Dataframe with stacked parameter df for each chain
# function getParamKernel(algo::MAlgoBGP)
function getParamCovariance(algo::MAlgoBGP)

	# select the last algo["past_iterations"] iterations from each chain
	lower_bound_index = maximum([1,algo.MChains[1].i-algo["past_iterations"]])

	# get all params from all chains 
	pardf = parameters(algo.MChains,lower_bound_index:algo.MChains[1].i)
	# |-------|----|------|----------|-----------|
	# | Row # | id | iter | a        | b         |
	# | 1     | 1  | 1    | 0.901685 | 0.983013  |
	# | 2     | 1  | 2    | 0.282522 | 0.763859  |
	# | 3     | 1  | 3    | 0.817773 | 0.268273  |
	# ...

	# compute Var-Cov matrix of parameters_to_sample
	# plus some small random noise
	VV = cov(array(pardf[:, algo.m.p2sample_sym])) + 0.0001 * Diagonal([1 for i=1:length(algo.m.p2sample_sym)])
	return VV

	# setup a MvNormal
	# MVN = MvNormal(VV)
	# return MVN
end


# function getNewCandidates!(algo::MAlgoBGP,MVN::MvNormal)
function getNewCandidates!(algo::MAlgoBGP,VV::Matrix)

	# update chain by chain
	for ch in 1:algo["N"]

		# TODO
		# getParamKernel could be in here
		# chain.temperature should parameterize MVN somehow (as in their toy example: multiplies the variance)
		# setup a MvNormal
		# VV2 = VV.*algo.MChains[ch].tempering
		VV2 = VV
		MVN = MvNormal(VV2)

		# shock parameters on chain index ch
		shock = rand(MVN) * algo.MChains[ch].shock_sd

		updateCandidateParam!(algo,ch,shock)

	end

end


function updateCandidateParam!(algo::MAlgoBGP,ch::Int,shock::Array{Float64,1})

	# get last param on that particular chain
	# as dataframe row
	oldpar = parameters(algo.MChains[ch],algo.MChains[ch].i-1)

	# add shock to each column of newpar
	newpar = copy(oldpar[algo.m.p2sample_sym])	# algo.m.p2sample_sym are in alphabetical order
	for c in 1:ncol(newpar)
		newpar[1,c] += shock[c]
	end

	# do bounds checking on newpar
	fitMirror!(newpar,algo.m.params_to_sample_df)

	# set as dict on algo.current_param
	fillinFields!(algo.candidate_param[ch],newpar)
end



