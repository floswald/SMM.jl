

# main dev routine for Mopt.jl

cd("/Users/florianoswald/git/MOpt.jl/")

include("src/mopt.jl")

# get a parameter vector
p = ["a" => 3.1 , "b" => 4.9]

# get some moments
mom = MomOpt.DataFrame(data=rand(3),sd=rand(3)*0.1,name=["alpha","beta","gamma"])

# which pars to sample
whichpar = MomOpt.DataFrame(name=["a","b"],lb=[-1,0],ub=[2,2])

# initiate
M = MomOpt.Mopt(p,whichpar,"Testobj",mom);