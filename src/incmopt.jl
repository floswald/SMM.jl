


type Mopt

	# chain setup
	iter             :: Int 	# iteration
	i                :: Int 	# ?
	use_last_run     :: Bool 	
	save_error       :: Bool
	params_to_sample :: Array{ASCIIString,1} 	# vector of param names to sample. i.e. fixed params don't appear here
	objfunc          :: ASCIIString 	# name of objective function
	run              :: Int 	# number of run
	shock_var        :: Float64 # variance of shock
	moments_to_use   :: Array{ASCIIString,1} 	#names of moments to use
	moments          :: DataFrame	# columns: data, sd
	np_shock         :: Float64
	save_freq        :: Int
	initial_value    :: Array{Float64,1} 	# initial parameter value
	params_all       :: Array{ASCIIString,1}	# names of all parameters
	param_descript   :: DataFrame 	# table with names and bounds of each paramtere
	N                :: Int 	# number of chains
	n_untempered     :: Int 	# number of chains untemprered

	# cluster setup
	mode         :: ASCIIString	# {'serial','mpi'}

	# paths for I/O
	file_chain       :: ASCIIString
	file_lastparam   :: ASCIIString
	file_errorparam  :: ASCIIString
	wd               :: ASCIIString
	source_on_nodes  :: ASCIIString
	logdir           :: ASCIIString


end