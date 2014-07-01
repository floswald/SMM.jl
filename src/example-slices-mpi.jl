

using MOpt


# initial value
p    = ["a" => 0.9 , "b" => -0.9]
# param bounds
pb   = [ "a" => [-1,1] , "b" => [-1,1] ]
# data moments
moms = DataFrame(moment=["alpha","beta"],data_value=[0.0,0.0],data_sd=rand(2))

# define a minization problem
mprob = MProb(p,pb,MOpt.objfunc_norm2,moms)

# compute slices in parallel
mprob_slices = MOpt.slices(mprob,30,true)
println(mprob_slices)



