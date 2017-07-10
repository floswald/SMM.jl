
module TestMProb

using Base.Test, MomentOpt, DataFrames, DataStructures




# Choose Algorithm and configure it

# malgo = MAlgoRandom(mprob,3)
# malgo["maxiter"]   = 100
# malgo["save_freq"] = 5
# malgo["mode"]      = "serial"


# TESTING MProb
# ==============

# test default constructor type
pb   = Dict( "a" => [0.5,0,1] , "b" => [0.6,0,1] )

moms = DataFrame(name=["alpha";"beta";"gamma"],value=[0.8;0.7;0.5],weight=rand(3))

@testset "Testing the MProb constructor" begin

	mprob = MProb()
	@test isa(mprob, MProb) == true

end


@testset "testing MProb methods" begin

	p   = OrderedDict(
		"a" => 0.1 , 
		"b" => 0.2 ,
		"c" => 0.5 )
	pb   = OrderedDict(
		"a" => [0.1; 0; 1] , 
		"b" => [0.2; 0; 1] )
	mprob = MProb();

	addParam!(mprob,p)	
	@test collect(MomentOpt.ps_names(mprob)) == Any[:a,:b,:c]
	@test length(mprob.params_to_sample) == 0

	mprob = MProb();
	addSampledParam!(mprob,pb)	
	@test collect(MomentOpt.ps_names(mprob)) == Any[:a,:b]
	@test collect(MomentOpt.ps2s_names(mprob)) == Any[:a,:b]

	mprob = MProb();
	addSampledParam!(mprob,"a",0.1,0,1)
	addSampledParam!(mprob,"b",0.1,0,1)
	addMoment!(mprob,moms)
	addEvalFunc!(mprob,MomentOpt.Testobj2)

	@test isa(mprob.objfunc,Function)

	@test collect(MomentOpt.ps_names(mprob)) == Any[:a,:b]
	@test collect(MomentOpt.ms_names(mprob)) == Any[:alpha,:beta,:gamma]

end



end # module 






