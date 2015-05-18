module TestBGPChain

using FactCheck, DataFrames, MOpt, Lazy



# TESTING Chains
# ==============
pb   = [ "a" => [3.1, 0,4] , "b" => [4.9,0,5,1] ]
moms = DataFrame(name=["alpha","beta","gamma"],value=[0.8,0.7,0.5],weight=rand(3))

facts("Testing BGPChain constructor") do
	
	mprob = @> MProb() addSampledParam!(pb) addMoment!(moms) addEvalFunc!(MOpt.Testobj)
	L = 9
	temp = 100.0
	shock = 12.0
	dist_tol = 0.001
	jumpprob = 0.05
	id = 180
	chain = BGPChain(id,mprob,L,temp,shock,dist_tol,jumpprob)

	@fact chain.i => 0 
	@fact chain.id => id

	context("length of members") do

		# test that all member except i are L long
		@fact length(chain.infos[:evals])  => L
		@fact length(chain.infos[:accept]) => L
		for nm in MOpt.ps2s_names(mprob)
			@fact nrow(chain.parameters) => L
		end
		for nm in MOpt.ms_names(mprob)
			@fact nrow(chain.moments) => L
		end
		@fact chain.tempering => temp
		@fact chain.shock_sd => shock
		@fact chain.jump_prob => jumpprob
		@fact chain.dist_tol => dist_tol
	end

	context("get correct dataframe labels") do
		@fact names(MOpt.parameters(chain)) => MOpt.ps2s_names(mprob)
		for nm in MOpt.ms_names(mprob)
			@fact nm in names(MOpt.moments(chain))[2:end] => true
		end
	end

end



end # module 