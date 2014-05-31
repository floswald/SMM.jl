

# test constructor

using Base.Test
using MomOpt


# get a parameter vector
p = ["a" => 3.1 , "b" => 4.9]

# get some moments
mom = DataFrame(data=rand(3),sd=rand(3)*0.1,name=["alpha","beta","gamma"])

# which pars to sample
whichpar = DataFrame(name=["a","b"],lb=[-1,0],ub=[2,2])

# initiate
M = Mopt(p,whichpar,"Testobj",mom);

@test isa(M,Mopt) == true


# test whether throws errors




# test dummy functions




