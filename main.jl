

# main dev routine for Mopt.jl

home = ENV["HOME"]
cd("$home/git/MOpt.jl")

# to develop with tests: run this
include("test/test_mopt.jl")

# to develop in main: run this
include("src/mopt.jl")





# workbench
# =========


# define a Test objective function
function Testobj(x::Dict,mom::DataFrame,whichmom::Array{ASCIIString,1})

	t0 = time()
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

	# status
	status = 1

	# time out
	t0 = time() - t0

	# return a dict
	ret = ["value" => v, "param" => x, "time" => t0, "status" => status, "moments" => mm]
	return ret
end


# get a parameter vector
p = ["a" => 3.1 , "b" => 4.9]
# define params to use with bounds
pb= [ "a" => [0,1] , "b" => [0,1] ]

# get some moments
moms = [
	"alpha" => [ 0.8 , 0.02 ],
	"beta"  => [ 0.8 , 0.02 ],
	"gamma" => [ 0.8 , 0.02 ]
]

submoms = ["alpha", "beta"]

# Define an Moment Optimization Problem
mprob = Mopt.MProb(p,pb,Testobj,moms;moments_subset=submoms)

# show
mprob


