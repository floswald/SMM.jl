	

using MOpt, Lazy

# data are generated from a bivariate normal
# with mu = [a,b] = [0,0]
# aim: 
# 1) sample [a',b'] from a space [-1,1] x [-1,1] and
# 2) find true [a,b] by computing distance(S([a',b']), S([a,b]))
#    and accepting/rejecting [a',b'] according to BGP
# 3) S([a,b]) returns a summary of features of the data

# initial value
pb    = ["a" => [1.9,-2,2] , "b" => [-0.9,-1,1] ] 
moms = DataFrame(name=["alpha","beta"],value=[0.0,0.0],weight=rand(2))
mprob = @> MProb() addSampledParam!(pb) addMoment(moms) addEvalFunc(MOpt.Testobj2)

# look at slice of the model: 
# how does the objective function behave 
# if we vary each parameter one by one, holding 
# the others fixed?

obj_slices = slices(mprob,30)

opts =[
	"N"               => 6,							# number of MCMC chains
	"maxiter"         => 500,						# max number of iterations
	"savefile"        => joinpath(pwd(),"MA.h5"),	# filename to save results
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


# run the estimation
runMOpt!(MA)


# save results
# save(MA,MA["savefile"])
	


