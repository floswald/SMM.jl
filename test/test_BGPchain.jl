module TestBGPChain

using FactCheck, DataFrames

include("../src/mopt.jl")


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
	
	mprob = Mopt.MProb(p,pb,Mopt.Testobj,moms)
	L = 9
	temp = 100.0
	shock = 12.0
	jumptol = 1.1
	id = 180
	chain = Mopt.BGPChain(id,mprob,L,temp,shock,jumptol)

	@fact chain.i => 0 
	@fact chain.id => id

	context("length of members") do

		# test that all member except i are L long
		@fact length(chain.infos["evals"])  => L
		@fact length(chain.infos["accept"]) => L
		for nm in Mopt.ps_names(mprob)
			@fact length(chain.parameters[nm]) => L
		end
		for nm in Mopt.ms_names(mprob)
			@fact length(chain.moments[nm]) => L
		end
		@fact chain.tempering => temp
		@fact chain.shock_sd => shock
		@fact chain.jumptol => jumptol
	end


end








end # module 