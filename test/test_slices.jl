module TestSlices

using FactCheck,DataFrames,MOpt,Lazy

	# initial value
	pb    = ["a" => [1.9,-2,2] , "b" => [-0.9,-1,1] ] 
	moms = DataFrame(name=["alpha","beta"],value=[0.0,0.0],weight=rand(2))
	mprob = @> MProb() addSampledParam!(pb) addMoment(moms) addEvalFunc(MOpt.Testobj2)

	# look at slice of the model: 
	# how does the objective function behave 
	# if we vary each parameter one by one, holding 
	# the others fixed?

	obj_slices = slices(mprob,30)


end # module 






