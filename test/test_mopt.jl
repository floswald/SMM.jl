home = ENV["HOME"]
cd("$home/git/MOpt.jl")
using Base.Test
using Mopt

include("src/mopt.jl")

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
p = [ "a" => 3.1 , 
			"b" => 4.9]

# define parameter bounds
pb= [ "a" => [0,1] , 
			"b" => [0,1] ]

# get some moments
moms = [
	"alpha" => [ 0.8 , 0.02 ],
	"beta"  => [ 0.8 , 0.02 ],
	"gamma" => [ 0.8 , 0.02 ]
]

# Define an Moment Optimization Problem
mprob = MProb(p,pb,moms,TestObj)

# Choose Algorithm and configure it

malgo = MAlgoRandom(mprob,3)
malgo["maxiter"]   = 100
malgo["save_freq"] = 5
malgo["mode"]      = "serial"









# TESTING
# =======

@test isa(mprob, MProb) == true

# test whether constructor throws errors
whichpar = Mopt.DataFrame(name=["a","c"],lb=[-1,0],ub=[2,2])
@test_throws ErrorException Moptim(p,whichpar,Testobj,mom; moments_to_use=["alpha"]);

whichpar = Mopt.DataFrame(name=["a","b"],lb=[-1,0],ub=[2,2])
mom = Mopt.DataFrame(data=rand(3),sd=rand(3)*0.1,NOTname=["alpha","beta","gamma"])
@test_throws ErrorException Moptim(p,whichpar,Testobj,mom; moments_to_use=["alpha"]);

# test dummy functions





