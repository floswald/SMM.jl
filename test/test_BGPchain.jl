module TestBGPChain

using FactCheck, MOpt



# TESTING Chains
# ==============
p    = ["a" => 3.1 , "b" => 4.9]
pb   = [ "a" => [0,1] , "b" => [0,1] ]
moms = [
	"alpha" => [ 0.8 , 0.02 ],
	"beta"  => [ 0.8 , 0.02 ],
	"gamma" => [ 0.8 , 0.02 ]
]



facts("Testing BGPChain constructor") do
	
	mprob = MProb(p,pb,Testobj,moms)
	L = 9
	temp = 100.0
	shock = 12.0
	accept_tol = 1.9
	dist_tol = 0.001
	jumpprob = 0.05
	id = 180
	chain = BGPChain(id,mprob,L,temp,shock,accept_tol,dist_tol,jumpprob)

	@fact chain.i => 0 
	@fact chain.id => id

	context("length of members") do

		# test that all member except i are L long
		@fact length(chain.infos[:evals])  => L
		@fact length(chain.infos[:accept]) => L
		for nm in MOpt.ps_names(mprob)
			@fact nrow(chain.parameters) => L
		end
		for nm in MOpt.ms_names(mprob)
			@fact nrow(chain.moments) => L
		end
		@fact chain.tempering => temp
		@fact chain.shock_sd => shock
		@fact chain.jump_prob => jumpprob
		@fact chain.dist_tol => dist_tol
		@fact chain.accept_tol => accept_tol
	end


end



end # module 