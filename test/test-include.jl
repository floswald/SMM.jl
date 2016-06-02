

using FactCheck,DataFrames,Lazy,MOpt
pb   = Dict("p1" => [0.0;-2;2] , "p2" => [0.0;-2;2] ) 
moms = DataFrame(name=["mu1";"mu2"],value=[0.0;0.0],weight=rand(2))
mprob_fail = @> MOpt.MProb() MOpt.addSampledParam!(pb) MOpt.addMoment!(moms) MOpt.addEvalFunc!(MOpt.Testobj_fails);




