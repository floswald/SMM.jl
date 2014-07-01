
# source on nodes
using MOpt


# initial value
p    = ["a" => 0.9 , "b" => -0.9]
# param bounds
pb   = [ "a" => [-1,1] , "b" => [-1,1] ]
# data moments
moms = DataFrame(moment=["alpha","beta"],data_value=[0.0,0.0],data_sd=rand(2))

# define a minization problem
mprob = MProb(p,pb,MOpt.objfunc_norm2,moms)

root = joinpath(ENV["HOME"],"git/MOpt.jl")

opts =[
	"N"               => 6,							# number of MCMC chains
	"mode"            => "mpi",						# mode: serial or mpi
	"maxiter"         => 500,						# max number of iterations
	"savefile"        => joinpath(root,"MA.h5"),	# filename to save results
	"source_on_nodes" => joinpath(root,"src/nodes.jl"),	
	"print_level"     => 1,							# increasing verbosity level of output
	"maxtemp"         => 100,						# tempering of hottest chain
	"min_shock_sd"    => 0.1,						# initial sd of shock on coldest chain
	"max_shock_sd"    => 1,							# initial sd of shock on hottest chain
	"past_iterations" => 30,						# num of periods used to compute Cov(p)
	"min_accept_tol"  => 100,						# ABC-MCMC cutoff for rejecting small improvements
	"max_accept_tol"  => 100,						# ABC-MCMC cutoff for rejecting small improvements
	"min_disttol"     => 0.1,						# distance tol for jumps from coldest chain
	"max_disttol"     => 0.1,						# distance tol for jumps from hottest chain
	"min_jump_prob"   => 0.05,						# prob of jumps from coldest chain
	"max_jump_prob"   => 0.2]						# prob of jumps from hottest chain

# setup the BGP algorithm
MA = MAlgoBGP(mprob,opts)