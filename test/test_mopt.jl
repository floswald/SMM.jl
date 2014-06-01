


using Base.Test
using Mopt



# test constructor

# get a parameter vector
p = ["a" => 3.1 , "b" => 4.9]

# get some moments
mom = Mopt.DataFrame(data=rand(3),sd=rand(3)*0.1,name=["alpha","beta","gamma"])

# which pars to sample
whichpar = Mopt.DataFrame(name=["a","b"],lb=[-1,0],ub=[2,2])

# initiate
M = Moptim(p,whichpar,"Testobj",mom);

@test isa(M,Moptim) == true


# test whether constructor throws errors
whichpar = Mopt.DataFrame(name=["a","c"],lb=[-1,0],ub=[2,2])
@test_throws ErrorException Moptim(p,whichpar,"Testobj",mom)


whichpar = Mopt.DataFrame(name=["a","b"],lb=[-1,0],ub=[2,2])
mom = Mopt.DataFrame(data=rand(3),sd=rand(3)*0.1,NOTname=["alpha","beta","gamma"])
@test_throws ErrorException Moptim(p,whichpar,"Testobj",mom)



# test dummy functions





