#'.. py:class:: MProb
#'
#'   the type MProb describes the optimization
#'   problem, it knows about the parameters
#'   to estimate, their names and bounds, 
#'   but also about the moments and their names.
type MProb

  # setup
  initial_value    :: Dict  # initial parameter value as a dict
  # params_to_sample :: Dict{ASCIIString,Array{Float64,1}}  # a dictionary of upper and lower bound for params we estimate (others are fixed)
  params_to_sample_df  :: DataFrame
  p2sample_sym     :: Array{Symbol,1} # column names of params to sample for dataframes
  objfunc          :: Function # objective function
  objfunc_opts     :: Dict     # options passed to the objective function, e.g. printlevel
  # moments          :: Dict  # a dictionary of moments to track
  moments          :: DataFrame # a dictionary of moments to track
  moments_subset   :: Array{ASCIIString}  # an array of moment names to subset objective funciton

  # constructor
  function MProb(
    initial_value,
    params_to_sample,
    objfunc,moments; 
    moments_subset=moments[:moment],objfunc_opts=Dict())

    this = new()

    # assert that moments has two entries for each moment: value and sd
    # @assert all(map(x -> length(x),values(moments)) .== 2)

    # assert that params_to_sample are subset of initial_value
    @assert issubset(keys(params_to_sample),keys(initial_value))

    # assert that params_to_sample has valid bounds: a vector with 2 increasing entries
    @assert all(map(x -> x[2]-x[1],values(params_to_sample)) .> 0)

    par2sample_name   = sort(collect(keys(params_to_sample)))
    par2sample_sym    = Symbol[x for x in par2sample_name]
    p0 = deepcopy(initial_value)

    pdf = DataFrame(param=collect(keys(params_to_sample)),lb=map(x -> x[1], values(params_to_sample)),ub=map(x -> x[2], values(params_to_sample)))
    sort!(pdf,cols=1)

    this.initial_value       = p0
    this.params_to_sample_df = pdf
    this.p2sample_sym        = par2sample_sym
    this.objfunc             = objfunc
    this.objfunc_opts        = objfunc_opts
    this.moments             = moments
    this.moments_subset      = moments_subset

    return(this)

  end # constructor

  # very simple constructor
  function MProb()
    this = new()
    this.initial_value       = Dict()
    this.params_to_sample_df = DataFrame()
    this.p2sample_sym        = []
    this.objfunc             = x -> x
    this.objfunc_opts        = Dict()
    this.moments             = DataFrame()
    this.moments_subset      = []
    return(this)
  end

end #type

function addParam!(m::MProb,name::ASCIIString,init)
  m.initial_value[name] = init
  return m 
end

function addSampledParam!(m::MProb,name::ASCIIString,init,lb,ub)
  m.initial_value[name] = init
  m.params_to_sample_df = rbind(m.params_to_sample_df,DataFrame(param=name,lb=lb,ub=ub))
  return m
end







#'.. py:function:: ps_names
#'
#'   returns the list of paramaters to sample
function ps_names(mprob::MProb)
  return(keys(mprob.initial_value))
end
function ps2s_names(mprob::MProb)
  return(mprob.p2sample_sym)
end

function ms_names(mprob::MProb)
  return(mprob.moments[:moment])
end

# evalute objective function
function evaluateObjective(m::MProb,p::Dict)

    x = eval(Expr(:call,m.objfunc,p,m.moments,m.moments_subset,m.objfunc_opts))
    gc()
    return x
end


#'.. py:function:: slices(m,pad)
#'
#'   computes slices for the objective function
function slices(m::MProb,npoints::Int,pad=0.1)

    # make a dict of grids for each param
    # loop over params!
    pdf = m.params_to_sample_df
    pranges = Dict{ASCIIString,Array{Float64,1}}()
    for irow in eachrow(pdf)
        lb = irow[:lb][1]
        ub = irow[:lb][1]
        pranges[irow[:param]] = linspace(irow[:lb][1], irow[:ub][1], npoints)  
    end
    # return pranges
    val_df = DataFrame()
    mom_df = DataFrame()
    for (k,v) in pranges
        println("currently computing slices over $k")
        dtmp = computeSlice(m,k,v)
        val_df = rbind(val_df,dtmp[1])
        mom_df = rbind(mom_df,dtmp[2])
    end
    return (val_df,mom_df)
   
end


#'.. py:function:: computeSlice(m,par,prange)
#'
#'   computes slices for the objective function
function computeSlice(m::MProb,par::ASCIIString,prange::Array{Float64,1})

    npar = length(prange)
    nmom = length(m.moments_subset)

    # make an array of different params
    # where par varies in prange
    pp = [deepcopy(m.initial_value) for i=1:npar]
    for i in 1:length(prange)
        pp[i][par] = prange[i]
    end

    v = pmap(x -> evaluateObjective(m,x), pp)

    mom_df = DataFrame(p_name = [par for i=1:(npar*nmom)], m_name=ASCIIString[ i for j=1:npar, i in m.moments_subset][:], p_val = repeat(prange,inner=[1],outer=[nmom]), m_val = zeros(npar*nmom))
    
    val_df = DataFrame(p_name = [par for i=1:(npar)], p_val = prange, f_val = zeros(npar), status = zeros(npar))

    for ip in 1:npar
        # fill in function values
        val_df[ip, :f_val ] = v[ip]["value"]
        val_df[ip, :status] = v[ip]["status"]

        for im in 1:nmom
            # fill in moments values
            mom_df[(mom_df[:p_val].==prange[ip]) & (mom_df[:m_name].==m.moments_subset[im]), :m_val ] = v[ip]["moments"][1,symbol(m.moments_subset[im])]
        end
    end

    return (val_df,mom_df)

end


function show(io::IO,m::MProb)


  print(io,"MProb Object:\n")
  print(io,"==============\n\n")
  print(io,"Parameters to sample:\n")
  print(io,m.params_to_sample_df)
  print(io,"\nMoment Table:\n")
  print(io,m.moments)
  print(io,"Moment to use:\n")
  print(io,m.moments_subset)
  print(io,"\n")
  print(io,"\nobjective function: $(m.objfunc)\n\n")
  # print(io,"\nobjective function call: $(Expr(:call,m.objfunc,m.current_param,m.moments,m.moments_to_use))\n")
  # if !m.prepared
  #   print(io,"\ncall MoptPrepare(m) to setup cluster\n")
  # else 
  #   print(io,"\nNumber of chains: $(m.N)\n")
  # end
  print(io,"\nEND SHOW MProb\n")
  print(io,"===========================\n")
end
