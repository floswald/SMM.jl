


type Mopt

	# chain setup
	initial_value    :: Dict 	# initial parameter value as a dict
	params_to_sample :: DataFrame 	
	objfunc          :: ASCIIString 	# name of objective function
	moments          :: DataFrame	# columns: data, sd
	moments_to_use   :: Array{ASCIIString,1} 	#names of moments to use
	shock_var        :: Float64 # variance of shock
	np_shock         :: Float64
	save_freq        :: Int
	# params_all       :: Array{ASCIIString,1}	# names of all parameters
	N                :: Int 	# number of chains
	n_untempered     :: Int 	# number of chains untemprered
	iter             :: Int 	# iteration
	run              :: Int 	# number of run
	i                :: Int 	# ?
	# use_last_run     :: Bool 	
	# save_error       :: Bool

	# cluster setup
	mode             :: ASCIIString	# {'serial','mpi'}

	# paths for I/O
	file_chain       :: ASCIIString
	file_lastparam   :: ASCIIString
	file_errorparam  :: ASCIIString
	wd               :: ASCIIString
	source_on_nodes  :: ASCIIString
	logdir           :: ASCIIString


	# constructor
	function Mopt(initial_value_,params_to_sample_,objfunc_,moments_; moments_to_use_=[""],shock_var_=1.0,np_shock_=1.0,save_freq_=25,N_=3,n_untempered_=1)

		iter = 0
		run = 0
		i   = 0

		# test format of moments_
		@assert haskey(moments_.colindex,:name)
		@assert haskey(moments_.colindex,:data)
		@assert haskey(moments_.colindex,:sd)

		# names in params_to_sample_ must be members of key(initial_value_)
		@assert length(setdiff(collect(keys(initial_value_)),params_to_sample_[:name])) == 0



		return new(initial_value_,params_to_sample_,objfunc_,moments_,moments_to_use_,shock_var_,np_shock_,save_freq_,N_,n_untempered_,iter,run,i)

	end	# constructor


end	#type

# tried that: doesnt work because the Mopt needs to take a generic parameter type, like a vector or a dict. otherwise every user will have their own param type.
# a test param type
# type Testparam
# 	a :: Float64
# 	b :: Float64

# 	function Testparam(a,b) 
# 	return new(a,b)
# 	end
# end


# define Test objective
function Testobj(x::Dict,mom::DataFrame,whichmom::Array{ASCIIString,1})

	mm = copy(mom)
	nm = names(mm)

	# add the model moments
	mm = cbind(mm,[x["a"] + x["b"] + rand() for i=1:nrow(mm)])
	names!(mm,[nm, :model])

	# subset to required moments only
	mm = mm[findin(mm[:name],whichmom),:]

	# compute distance
	v = sum((mm[:data]-mm[:model])^2)
end
