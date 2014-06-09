

# main dev routine for Mopt.jl

home = ENV["HOME"]
cd("$home/git/MOpt.jl")

# to develop with tests: run this
include("test/test_mopt.jl")

# to develop in main: run this
include("src/mopt.jl")





# workbench
# =========





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

# call objective
x = Mopt.Testobj(p,moms,submoms)

# Define an Moment Optimization Problem
mprob = Mopt.MProb(p,pb,Testobj,moms;moments_subset=submoms)

# show
mprob


