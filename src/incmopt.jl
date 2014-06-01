


type Moptim

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
	paths :: Dict

	prepared :: Bool

	# constructor
	function Moptim(initial_value,params_to_sample,objfunc,moments; moments_to_use=moments[:name],shock_var=1.0,np_shock=1.0,save_freq=25,N=3,n_untempered=1,mode="serial",paths=["chain"=> ".","lastparam"=>".","errorparam"=>".","wd"=>".","source_on_nodes"=>"workers.jl","logdir"=>"./log"])

		iter = 0
		run = 0
		i   = 0

		# test format of moments_
		@assert haskey(moments.colindex,:name)
		@assert haskey(moments.colindex,:data)
		@assert haskey(moments.colindex,:sd)

		# names in params_to_sample_ must be members of key(initial_value_)
		@assert length(setdiff(collect(keys(initial_value)),params_to_sample[:name])) == 0

		if mode == "serial"
			prepared = true
		else
			prepared = false
		end


		return new(initial_value,params_to_sample,objfunc,moments,moments_to_use,shock_var,np_shock,save_freq,N,n_untempered,iter,run,i,mode,paths,prepared)

	end	# constructor
end	#type

function show(io::IO,m::Moptim)
	print(io,"Moptim Object:\n")
	print(io,"==============\n\n")
	print(io,"Parameters to sample:\n")
	print(io,m.params_to_sample)
	print(io,"\nMoment Table:\n")
	print(io,m.moments)
	print(io,"\nMoment to use:\n")
	print(io,m.moments_to_use)
	print(io,"\n\nMode: $(m.mode)\n")
	if !m.prepared
		print(io,"\n\ncall MoptPrepare(m) to setup cluster\n")
	else 
		print(io,"Number of chains: $(m.N)\n")
	end
end


# tried that: doesnt work because the Mopt needs to take a generic parameter type, like a vector or a dict. otherwise every user will have their own param type. I don't even know the name of that user type, let alone it's structure.
# a test param type
# type Testparam
# 	a :: Float64
# 	b :: Float64

# 	function Testparam(a,b) 
# 	return new(a,b)
# 	end
# end



# call this function to set up
# the cluster. 
# maps choice of "mode" into an action
function MoptPrepare(m::Moptim)

	if m.prepared
		println("you're good to go")
		return nothing
	end

	# setup cluster

	# set m.N = # of workers

	# send data to workers

	# set m.prepared = true

end





# define Test objective function
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
