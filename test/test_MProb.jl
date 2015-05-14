
module TestMProb

using FactCheck, MOpt




# Choose Algorithm and configure it

# malgo = MAlgoRandom(mprob,3)
# malgo["maxiter"]   = 100
# malgo["save_freq"] = 5
# malgo["mode"]      = "serial"


# TESTING MProb
# ==============

# test default constructor type
pb   = [ "a" => [0.5,0,1] , "b" => [0.6,0,1] ]

moms = DataFrame(name=["alpha","beta","gamma"],value=[0.8,0.7,0.5],weight=rand(3))

facts("Testing the MProb constructor") do

	mprob = MProb()
	@fact isa(mprob, MProb) => true

end


facts("testing MProb methods") do

	pb   = [ 
		"a" => [0.1, 0, 1] , 
		"b" => [0.2, 0, 1] ]
	mprob = MProb();
	addSampledParam!(mprob,pb)	
	@fact sort(collect(MOpt.ps_names(mprob))) == sort({:a,:b}) => true

	mprob = MProb();
	addSampledParam!(mprob,"a",0.1,0,1)
	addSampledParam!(mprob,"b",0.1,0,1)
	addMoment!(mprob,moms)
	addEvalFunc!(mprob,MOpt.Testobj2)

	@fact isa(mprob.objfunc,Function) => true

	@fact sort(collect(MOpt.ps_names(mprob))) == sort({:a,:b}) => true
	@fact sort(collect(MOpt.ms_names(mprob))) == sort({:alpha,:beta,:gamma}) => true

end


end # module 






