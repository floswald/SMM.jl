# the type MProb describes the optimization
# problem, it knows about the parameters
# to estimate, their names and bounds, 
# but also about the moments and their names.
type MProb

  # setup
  initial_value    :: Dict  # initial parameter value as a dict
  # params_to_sample :: Dict{ASCIIString,Array{Float64,1}}  # a dictionary of upper and lower bound for params we estimate (others are fixed)
  params_to_sample_df  :: DataFrame
  p2sample_sym     :: Array{Symbol,1} # column names of params to sample for dataframes
  objfunc          :: Function # objective function
  # moments          :: Dict  # a dictionary of moments to track
  moments          :: DataFrame # a dictionary of moments to track
  moments_subset   :: Array{ASCIIString}  # an array of moment names to subset objective funciton

  # constructor
  function MProb(
    initial_value,
    params_to_sample,
    objfunc,moments; 
    moments_subset=moments[:moment])

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

    return new(p0,pdf,par2sample_sym,objfunc,moments,moments_subset)

  end # constructor
end #type

# returns the list of paramaters to sample
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

    x = eval(Expr(:call,m.objfunc,p,m.moments,m.moments_subset))
    return x
end


function computeSlice(m::MProb,npoints::Int,par::ASCIIString,pad)

    pdf = m.params_to_sample_df
    # make a deepcopy of initial_value of parameters
    pp = deepcopy(m.initial_value)
    # generate a param range
    lb = pdf[pdf[:param].==par,:lb][1]
    ub = pdf[pdf[:param].==par,:ub][1]

    prange = linspace( lb+(ub-lb)*pad/2, lb+(ub-lb)*(1-pad/2),npoints)
    # make first row for this par
    pp[par] = prange[1]
    v = evaluateObjective(m,pp)
    df = DataFrame(p_name=par,p_val=prange[1],f_val=v["value"])

    for i in 2:npoints
        pp[par] = prange[i]
        v = evaluateObjective(m,pp)
        push!(df,{par,prange[i],v["value"]})
    end
    return df
end

function slices(m::MProb,npoints::Int,parallel::Bool,pad=0.1)

    npar = nrow(m.params_to_sample_df)
    if parallel
        v = pmap( x -> computeSlice(m,npoints,x,pad), m.params_to_sample_df[:param])
     else
        v =  map( x -> computeSlice(m,npoints,x,pad), m.params_to_sample_df[:param])
    end
    
    df = v[1]
    if npar>1
        for i in 2:npar
            df = rbind(df,v[i])    
        end
    end

    return df
end


# function slices2(m::MProb,npoints::Int,parallel::Bool,pad=0.1)

#     # make a dict of grids for each param
#     pdf = m.params_to_sample_df
#     pranges = Dict{ASCIIString,Array{Float64,1}}[]
#     for irow in eachrow(pdf)
#         lb = irow[:lb][1]
#         ub = irow[:lb][1]
#         push!(pranges,[irow[:param] => linspace(irow[:lb][1], irow[:ub][1], npoints) ] )
#     end
#     dd = computeSlice2()

   
# end

# function computeSlice2(m::MProb,prange::Dict{ASCIIString,Array{Float64,1}},par::ASCIIString,pad)

#     pdf = m.params_to_sample_df
#     # make a deepcopy of initial_value of parameters
#     pp = deepcopy(m.initial_value)

#     # make first row for this par
#     pp[par] = prange[1]
#     v = evaluateObjective(m,pp)
#     df = DataFrame(p_name=par,p_val=prange[1],f_val=v["value"])

#     for i in 2:npoints
#         pp[par] = prange[i]
#         v = evaluateObjective(m,pp)
#         push!(df,{par,prange[i],v["value"]})
#     end
#     return df
# end




function plotSlices(m::MProb,x::DataFrame)
    npars = length(m.p2sample_sym)
    nrows = floor(sqrt(npars))
    ncols = ceil(npars/nrows)
    pid = 0
    for sdf in groupby(x, :p_name)
        subplot(nrows,ncols,pid)
        pid += 1
        plot(sdf[:p_val],sdf[:f_val])
        title(sdf[1,:p_name])
    end
    suptitle("slices")
end



# gadfly works, but gets mixed up with PyPlot.plot!
# function plotSlices(x::DataFrame,filename)

#     p=plot(x,xgroup="p_name",x="p_val",y="f_val",Geom.subplot_grid(Geom.line))
#     draw(PDF(filename),p)

# end








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