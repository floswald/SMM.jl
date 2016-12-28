


module TestAlgoBGP

using Base.Test, DataFrames, MOpt

pb   = Dict( "a" => [0.3; 0;1] , "b" => [0.4;0;1])
moms = DataFrame(name=["alpha","beta","gamma"],value=[0.8;0.7;0.5],weight=rand(3))
mprob = MProb() 
addSampledParam!(mprob,pb) 
addMoment!(mprob,moms) 
addEvalFunc!(mprob,MOpt.Testobj)

@testset "AlgoBGP" begin

	@testset "Constructor" begin

		MA = MAlgoBGP(mprob)

		@test isa(MA,MAlgo) == true
		@test isa(MA,MAlgoBGP) == true
		@test isa(MA.m,MProb) == true

		@test MA.i == 0
		@test length(MA.chains) == 3

		for ix = 1:length(MA.chains)
			@test MA.current_param[ix] == mprob.initial_value
		end

	end

end



# @testset "testing getNewCandidates(MAlgoBGP)" begin

# 	@testset "test fitMirror(x,ub,lb)" begin

# 		x = 12.1
# 		ub = 12.0
# 		lb = -1.1
# 		z = MOpt.fitMirror(x,lb,ub)
# 		@test (z < ub) & (z>lb) == true

# 		x = rand()
# 		z = MOpt.fitMirror(x,lb,ub)
# 		@test z-x<eps() == true

# 		x = -1.11
# 		z = MOpt.fitMirror(x,lb,ub)
# 		@test (z < ub) & (z>lb) == true

# 	end

# 	@testset "test fitMirror!(df,dict)" begin

# 		df = DataFrame(a=1.2,b=-1.5)
# 		MOpt.fitMirror!(df,mprob.params_to_sample)
# 		@test df[1][1] > pb["a"][2] == true
# 		@test df[1][1] < pb["a"][3] == true
# 		@test df[2][1] > pb["b"][2] == true
# 		@test df[2][1] < pb["b"][3] == true

# 	end

# 	MA = MAlgoBGP(mprob,opts)

# 	# fill chains with random values
# 	for iter =1:(MA["maxiter"]-1)
# 		# update all chains to index 1
# 		MOpt.incrementChainIter!(MA.MChains)

# 		# set a parameter vector on all chains using appendEval!
# 		myp = Dict(:a => rand() , :b => rand())
# 		mym = Dict( :alpha => rand(), :beta => rand(), :gamma => rand() )
# 		ev = Eval(); ev.value   =  1.1; ev.params  = myp; ev.time = 0; ev.status  = 1
# 		setMoment(ev,mym)
# 		for ich in 1:MA["N"]
# 			MOpt.appendEval!(MA.MChains[ich],ev,true,1.0)
# 		end
# 	end

# 	@testset "test getParamCovariance" begin  
# 		# get all parameters
# 		lower_bound_index = maximum([1,MA.MChains[1].i-MA["past_iterations"]])
# 		pars = parameters(MA.MChains,lower_bound_index:MA.MChains[1].i)

# 		# get the last MA["past_iterations"] iterations from each chain
# 		# get parameter_to_sample names as symbols 

# 		# get covariance matrix of those
# 		VV = cov( convert(Array,pars[:, ]) ) + 0.0001 * Diagonal([1 for p in pb])

# 		# get kernel and check VV
# 		myVV = MOpt.getParamCovariance(MA)

# 		@test myVV[:] == VV[:] == true

# 	end

# 	@testset "test updateCandidateParam!" begin

# 		# get a "kernel"/Covariance matrix
# 		# pos def matrix:
# 		myVV = rand(length( MOpt.ps2s_names(MA.m)),length(MOpt.ps2s_names(MA.m)))
# 		myVV = myVV * myVV'

# 		# manually create next parameter guess on
# 		# each chain
# 		newp = Array(DataFrame,MA["N"])

# 		# on all chains, the parameter entry number ix+1 
# 		# must correspond to entry number ix plus a shock from the kernel
# 		ich = 1; ix = 10
# 		MVN = MOpt.MvNormal( myVV.*MA.MChains[ich].tempering )
# 		shockd = Dict(zip(MOpt.ps2s_names(MA.m),rand(MVN) ))

# 		eval_before = getLastEval(MA.MChains[ich]).params # get a dataframe of row ix-1
# 		MOpt.jumpParams!(MA,ich,shockd)

# 		for (k,v) in eval_before
# 			@test MA.current_param[ich][k] == fitMirror(v + shockd[k], MA.m.params_to_sample[k])
# 		end
# 	end
# end


# @testset "testing swaprows" begin

# 	@testset "testing swapRows!" begin

# 		MA = MAlgoBGP(mprob,opts)

# 		# fill chains with random values up to iteration ix
# 		ix = 5
# 		for iter =1:ix
# 			# update all chains to index 1
# 			MOpt.incrementChainIter!(MA.MChains)

# 			# set a parameter vector on all chains using appendEval!
# 			for ich in 1:MA["N"]
# 				myp = Dict(:a => rand() , :b => rand())
# 				mym = Dict( :alpha => rand(), :beta => rand(), :gamma => rand() )
# 				ev = Eval(); ev.value   =  1.1; ev.params  = myp; ev.time = 0; ev.status  = 1
# 				setMoment(ev,mym)
# 				MOpt.appendEval!(MA.MChains[ich],ev,true,1.0)
# 			end

# 		end

# 		# get a pair (i,j) of chains
# 		pair = Pair(1,3)

# 		# get params and moms
# 		p1 = parameters(MA.MChains[pair.first ],ix)
# 		p2 = parameters(MA.MChains[pair.second],ix)
# 		m1 = MA.MChains[pair.first ].moments[ix,MA.MChains[pair.first ].moments_nms]
# 		m2 = MA.MChains[pair.second].moments[ix,MA.MChains[pair.second].moments_nms]
# 		v1 = evals(MA.MChains[pair.first ],ix)
# 		v2 = evals(MA.MChains[pair.second],ix)


# 		# exchange
# 		MOpt.swapRows!(MA,pair,ix)

# 		# check parameters
# 		@test parameters(MA.MChains[pair.first ],ix) == p2 == true
# 		@test parameters(MA.MChains[pair.second],ix) == p1 == true
# 		@test MA.MChains[pair.first ].moments[ix,MA.MChains[pair.first ].moments_nms] == m2 == true
# 		@test MA.MChains[pair.second].moments[ix,MA.MChains[pair.second].moments_nms] == m1 == true
# 		@test evals(MA.MChains[pair.first ],ix) == v2 == true
# 		@test evals(MA.MChains[pair.second],ix) == v1 == true

# 	end

	
# end

# @testset "testing accept reject" begin


# 	@testset "testing initial period" begin
# 		MA = MAlgoBGP(mprob,opts)

# 		# get a return value
# 		v = map( x -> MOpt.evaluateObjective(MA.m,x), MA.current_param)

# 		# set iteration on chains and algo = 1
# 		MA.i = 1
# 		MOpt.incrementChainIter!(MA.MChains)
# 		MOpt.doAcceptRecject!(MA,v)

# 		# all accepted: 
# 		@test all(infos(MA.MChains,1)[:accept]) == true
# 		@test all(infos(MA.MChains,1)[:status] .== 1.0) == true

# 		# all params equal to initial value
# 		for nm in MOpt.ps2s_names(MA)
# 			#print( " $nm $(parameters(MA,1,1,nm)) $(pb[string(nm)][1] ) ")
# 			@test parameters(MA,1,1,nm) == roughly(pb[string(nm)][1])
# 		end

# 		# next iteration
# 		# set options
# 		MA.i = 2
# 		MOpt.incrementChainIter!(MA.MChains)
# 		MOpt.doAcceptRecject!(MA,v)
# 	end

# 	@testset "testing whether params get accept/rejected" begin

# 		# first iteration
# 		MA = MAlgoBGP(mprob,opts)
# 		MA.i = 1
# 		MOpt.incrementChainIter!(MA.MChains)
# 		v = map( x -> MOpt.evaluateObjective(MA.m,x), MA.current_param)

# 		MOpt.doAcceptRecject!(MA,v)

# 		# second iteration
# 		MA.i = 2
# 		MOpt.incrementChainIter!(MA.MChains)

# 		# get randome parameters
# 		for i in 1:length(v) 
# 			for (k,va) in MA.current_param[i]
# 				MA.current_param[i][k] = va + rand()
# 			end
# 		end

# 		# get 2 return values
# 		v1 = map( x -> MOpt.evaluateObjective(MA.m,x), MA.current_param)
# 		v2 = map( x -> MOpt.evaluateObjective(MA.m,x), MA.current_param)

# 		# assign a very good bad value on each chain: 
# 		for i in 1:length(v) 
# 			v1[i].value = v1[i].value - 100.0
# 			setMoment(v1[i],  Dict( :alpha => rand(), :beta => rand(), :gamma => rand() ))
# 			v2[i].value = v2[i].value + 100.0
# 			setMoment(v2[i],  Dict( :alpha => rand(), :beta => rand(), :gamma => rand() ) )
# 		end
# 		MAs = MAlgo[]
# 		push!(MAs,deepcopy(MA),deepcopy(MA))
# 		MOpt.doAcceptRecject!(MAs[1],v1)
# 		MOpt.doAcceptRecject!(MAs[1],v2)

# 		for ma in MAs
# 			for ch in 1:ma["N"]
# 				ev= getEval(ma.MChains[ch],ma.i)
# 				Last_ev= getEval(ma.MChains[ch],ma.i-1)

# 				if infos(ma.MChains[ch],ma.i)[:accept][1]
# 					# if accepted, the parameter vector in ev is equal to current_param on that chain
# 					ev_p = paramd(ev)
# 					ev_m = ev.simMoments
# 					for (k,v) in ev_p
# 						@test v == ma.current_param[ch][k]
# 					end
# 				else
# 					# if not, parameters(ma.MChains[ch],ma.i) == parameters(ma.MChains[ch],ma.i-1)
# 					ev_p = paramd(Last_ev)
# 					ev_m = Last_ev.simMoments
# 					for (k,v) in ev_p
# 						@test v == not(ma.current_param[ch][k])
# 					end
# 				end
# 			end
# 		end
# 	end
# end

# @testset "testing exchangeMoves" begin

# 	@testset "testing exchange Moves: some chains exchange" begin

# 	println("this test sets min_jump_prob=1 to make sure that there is a jump for sure")

# 		newopts =Dict("N"=>20,"shock_var"=>1.0,"mode"=>"serial","maxiter"=>100,"path"=>".","maxtemp"=>100,"min_shock_sd"=>0.1,"max_shock_sd"=>1.0,"past_iterations"=>30,"min_disttol"=>0.01,"max_disttol"=>0.01,"min_jump_prob"=>1.0,"max_jump_prob"=>1.0,"print_level" => 3) 

# 		MA = MAlgoBGP(mprob,newopts)
# 		MA.i = 1

# 		MOpt.incrementChainIter!(MA.MChains)
# 		v = map( x -> MOpt.evaluateObjective(MA.m,x), MA.current_param)
# 		for ch in 1:MA["N"] MOpt.appendEval!(MA.MChains[ch],v[ch],true,1.0) end

# 		MA.i = 2
# 		MOpt.incrementChainIter!(MA.MChains)

# 		# we set chain 2 very close to chain 
# 		# to make a jump
# 		ch1 = 1
# 		ch2 = 2
# 		v[ch2].value = v[ch1].value
# 		MOpt.appendEval!(MA.MChains[ch1],v[ch1],false,1.0) 
# 		MOpt.appendEval!(MA.MChains[ch2],v[ch2],false,1.0) 

# 		for i in 2:length(v) 
# 			# for all other chains set jump_prob = 0 
# 			MA.MChains[i].jump_prob = 0.0
# 			MOpt.appendEval!(MA.MChains[i],v[i],false,1.0) 
# 		end

# 		# I check that the current value for both chains are very close
# 		v1 = evals(MA.MChains[ch1],MA.MChains[ch1].i)[1]
# 		v2 = evals(MA.MChains[ch2],MA.MChains[ch2].i)[1]

# 		MOpt.exchangeMoves!(MA)

# 		# expect we exchanged ch1 and ch2 but nothing else
# 		ex_1 = infos(MA.MChains[ch1],MA.i)[:exchanged_with][1]
# 		@test  ex_1 != 0 == true
# 		for i in 1:length(v) 
# 			if ((i != ch1)  & ( i != ex_1) )
# 				@test infos(MA.MChains[i],MA.i)[:exchanged_with][1] == 0 ==true
# 			else
# 				# dont' know the pairs!
# 			end
# 		end

# 	end

# end


# @testset "testing reading and saving of algo" begin

# 	MA = MAlgoBGP(mprob,opts)
# 	runMOpt!(MA)
# 	save(MA,joinpath(pwd(),"test.h5"))

#     dd = Dict()
# 	for (k,v) in MA.opts
# 		dd[k] = "$v"
# 	end

# 	x = MOpt.readAlgoBGP(joinpath(pwd(),"test.h5"))
# 	for (k,v) in dd
# 		@test haskey(dd,k) == true
# 		@test v == dd[k] == true
# 	end

# 	ich = rand(1:MA["N"]) 	# pick a random chain
# 	# cc = @where(x["parameters"],:chain_id .== ich)
# 	@test all(parameters(MA.MChains[ich])[:a] .- x["params"][x["params"][:chain_id] .== ich, :a] .< 1e-6) == true
# 	@test all(parameters(MA.MChains[ich])[:b] .- x["params"][x["params"][:chain_id] .== ich, :b] .< 1e-6) == true

# 	rm(joinpath(pwd(),"test.h5"))

# end


# @testset "testing intermittent saving of algo" begin

# 	opts["filename"] = joinpath(pwd(),"test.h5")
# 	opts["save_frequency"] = 2
# 	MA = MAlgoBGP(mprob,opts)

# 	# replicate run()
# 	for i in 1:MA["maxiter"]

# 		MA.i = i
# 		MOpt.computeNextIteration!( MA )

# 		if mod(i,MA["save_frequency"]) == 0
# 			save(MA,MA["filename"])
# 			# the saved dataset must now contain 
# 			# all data up to iteration i
# 			ff5 = MOpt.h5open(MA["filename"],"r")
# 			ich = rand(1:MA["N"]) 	# pick a random chain
# 			chain_a = read(ff5,"chain/$ich/parameters/a")
# 			chain_b = read(ff5,"chain/$ich/parameters/b")
# 			close(ff5)

# 			@test length(find(chain_a.!=0)) == i
# 			@test length(find(chain_b.!=0)) == i

# 			@test all(parameters(MA.MChains[ich])[:a] .- chain_a .< 1e-6) == true
# 			@test all(parameters(MA.MChains[ich])[:b] .- chain_b .< 1e-6) == true
# 		end
# 	end

# 	rm(MA["filename"])

# end





end # module

