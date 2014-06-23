


# defining an abstract chain type in case 
# some algorithm need additional information
# on top of (eval, moment, parameter)
abstract AbstractChain


# methods for a single chain
# ==========================

parameters(c::AbstractChain, i::UnitRange{Int}) = c.parameters[i,:]
parameters(c::AbstractChain, i::Int)            = parameters(c, i:i)
parameters(c::AbstractChain)                    = c.parameters
moments(c::AbstractChain)                       = c.moments[i,:]
moments(c::AbstractChain, i::UnitRange{Int})    = c.moments[i,:]
moments(c::AbstractChain, i::Int)               = moments(c, i:i)
infos(c::AbstractChain)                         = c.infos
infos(c::AbstractChain, i::UnitRange{Int})      = c.infos[i,:]
infos(c::AbstractChain, i::Int)                 = infos(c, i:i)
evals(c::AbstractChain, i::UnitRange{Int})      = c.infos[i,:evals]
evals(c::AbstractChain, i::Int)                 = evals(c, i:i)
allstats(c::AbstractChain,i::UnitRange{Int})        = cbind(c.infos[i,:],c.parameters[i,:],c.moments[i,:])
allstats(c::AbstractChain,i::Int)                   = cbind(infos(c, i:i),parameters(c, i:i),moments(c, i:i))

# appends values from objective function
# at CURRENT iteration
function appendEval!(chain::AbstractChain, vals::Dict, ACC::Bool, status::Int, prob::Float64)
    chain.infos[chain.i,:evals] = vals["value"]
    chain.infos[chain.i,:prob] = prob
    chain.infos[chain.i,:accept] = ACC
    chain.infos[chain.i,:status] = status
    for im in chain.moments_nms
        chain.moments[chain.i,im] = vals["moments"][string(im)][1]
    end
    for ip in chain.params_nms
        chain.parameters[chain.i,ip] = vals["params"][string(ip)][1]
    end
  return nothing
end

# methods for an array of chains
# ==============================

# TODO ideally dispatch on MC::Array{AbstractChain}
# but that doesn't work. errors with
# no method for parameters(BGPChain)
#
# return an rbind of params from all chains
function parameters(MC::Array,i::UnitRange{Int})
    if !isa(MC[1],AbstractChain)
        error("must give array of AbstractChain") 
    end
    r = parameters(MC[1],i) 
    if length(MC)>1
        for ix=2:length(MC)
            r = rbind(r,parameters(MC[ix],i))
        end
    end
    return r
end

function parameters(MC::Array,i::Int)
    parameters(MC,i:i)
end
function parameters(MC::Array)
    parameters(MC,1:MC[1].i)
end

function moments(MC::Array,i::UnitRange{Int})
    if !isa(MC[1],AbstractChain)
        error("must give array of AbstractChain") 
    end
    r = moments(MC[1],i) 
    if length(MC)>1
        for ix=2:length(MC)
            r = rbind(r,moments(MC[ix],i))
        end
    end
    return r
end

function moments(MC::Array)
    moments(MC,1:MC[1].i)
end

function infos(MC::Array,i::Int)
    infos(MC,i:i)
end

function infos(MC::Array,i::UnitRange{Int})
    if !isa(MC[1],AbstractChain)
        error("must give array of AbstractChain") 
    end
    r = infos(MC[1],i) 
    if length(MC)>1
        for ix=2:length(MC)
            r = rbind(r,infos(MC[ix],i))
        end
    end
    return r
end

function infos(MC::Array)
    infos(MC,1:MC[1].i)
end

function evals(MC::Array,i::UnitRange{Int})
    if !isa(MC[1],AbstractChain)
        error("must give array of AbstractChain") 
    end
    r = infos(MC[1],i) 
    if length(MC)>1
        for ix=2:length(MC)
            r = rbind(r,evals(MC[ix],i))
        end
    end
    return r
end

function evals(MC::Array)
    evals(MC,1:MC[1].i)
end



function allstats(MC::Array) 
    cbind(infos(MC),parameters(MC)[MC[1].params_nms],moments(MC)[MC[1].moments_nms])
end



# the default chain type
# we create a dictionary with arrays
# for each parameters
type Chain <: AbstractChain
  i::Int              # current index
  infos      ::DataFrame   
  parameters ::DataFrame  
  moments    ::DataFrame 
  params_nms ::Array{Symbol,1}  # DataFrame names of parameters (i.e. exclusive of "id" or "iter", etc)
  moments_nms::Array{Symbol,1}  # DataFrame names of moments

  function Chain(MProb,L)
    infos      = DataFrame(iter=1:L, evals = zeros(Float64,L), accept = zeros(Bool,L), status = zeros(Int,L), exhanged_with=zeros(Int,L), prob=zeros(Float64,L))
    par_nms = Symbol[ symbol(x) for x in ps_names(MProb) ]
    mom_nms = Symbol[ symbol(x) for x in ms_names(MProb) ]
    parameters = convert(DataFrame,zeros(L,length(par_nms)+1))
    moments    = convert(DataFrame,zeros(L,length(mom_nms)+1))
    names!(parameters,[:iter, par_nms])
    names!(moments   ,[:iter, mom_nms])
    return new(0,infos,parameters,moments,par_nms,mom_nms)
  end
end


# update the iteration count on each chain
function updateIterChain!(MC::Array)
    for ix in 1:length(MC)
        MC[ix].i += 1
    end 
end

# saves a chain to a HF5 file at a given path
function saveToHDF5(chain::AbstractChain, ff5::HDF5File, path::ASCIIString)
    simpleDataFrameSave(chain.parameters,ff5, "$path/parameters")
    simpleDataFrameSave(chain.infos,ff5, "$path/infos")
    simpleDataFrameSave(chain.moments,ff5, "$path/moments")
end

function simpleDataFrameSave(dd::DataFrame,ff5::HDF5File, path)
    for nn in names(dd)
        write(ff5,"$path/$(string(nn))",convert(Array{Float64,1},dd[nn]))
    end
end


