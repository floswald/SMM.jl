


using Base.Test
using Mopt



# test constructor# define a Test objective function
function Testobj(x::Dict,mom::DataFrame,whichmom::Array{ASCIIString,1})

	mm = copy(mom)
	nm = names(mm)

	# perform computations of the objective function
	# and add the resulting model moments
	mm = cbind(mm,[x["a"] + x["b"] + i for i=1:nrow(mm)])
	names!(mm,[nm, :model])

	# subset to required moments only
	mm = mm[findin(mm[:name],whichmom),:]

	# compute distance
	v = sum((mm[:data].-mm[:model]).^2)
end

# get a parameter vector
p = ["a" => 3.1 , "b" => 4.9]

# get some moments
mom = DataFrame(data=rand(3),sd=rand(3)*0.1,name=["alpha","beta","gamma"])

# which pars to sample
whichpar = Mopt.DataFrame(name=["a","b"],lb=[-1,0],ub=[2,2])

# initiate
M = Moptim(p,whichpar,Testobj,mom; moments_to_use=["alpha"]);

@test isa(M,Moptim) == true


# test whether constructor throws errors
whichpar = Mopt.DataFrame(name=["a","c"],lb=[-1,0],ub=[2,2])
@test_throws ErrorException Moptim(p,whichpar,Testobj,mom; moments_to_use=["alpha"]);


whichpar = Mopt.DataFrame(name=["a","b"],lb=[-1,0],ub=[2,2])
mom = Mopt.DataFrame(data=rand(3),sd=rand(3)*0.1,NOTname=["alpha","beta","gamma"])
@test_throws ErrorException Moptim(p,whichpar,Testobj,mom; moments_to_use=["alpha"]);



# test dummy functions





