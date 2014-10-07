


# defining an abstract chain type in case 
# some algorithm need additional information
# on top of (eval, moment, parameter)
abstract AbstractChain


# the default chain type
# we create a dictionary with arrays
# for each parameters
type Chain <: AbstractChain
    i            :: Int              # current index
    infos        :: DataFrame
    parameters   :: DataFrame
    moments      :: DataFrame
    params_nms   :: Array{Symbol,1}  # DataFrame names of parameters (i.e. exclusive of "id" or "iter", etc)
    params2s_nms :: Array{Symbol,1}  # DataFrame names of parameters to sample
    moments_nms  :: Array{Symbol,1}  # DataFrame names of moments
  
    function Chain(MProb,L)
        # infos      = DataFrame(iter=0, evals = 0.0, accept = true, status = 0, exhanged_with=0, prob=0.0)
        infos      = DataFrame(iter=1:L, evals =zeros(Float64,L), accept = zeros(Bool,L), status = zeros(Int,L), exhanged_with=zeros(Int,L), prob=zeros(Float64,L), eval_time=zeros(Float64,L))
        par_nms    = Symbol[ symbol(x) for x in ps_names(MProb) ]
        par2s_nms  = Symbol[ symbol(x) for x in ps2s_names(MProb) ]
        mom_nms    = Symbol[ symbol(x) for x in ms_names(MProb) ]
        parameters = convert(DataFrame,zeros(L,length(par_nms)+1))
        moments    = convert(DataFrame,zeros(L,length(mom_nms)+1))
        names!(parameters,[:iter, par_nms])
        names!(moments   ,[:iter, mom_nms])
        return new(0,infos,parameters,moments,par_nms,par2s_nms,mom_nms)
    end
end

# methods for a single chain
# ==========================

function parameters(c::AbstractChain, i::Union(Integer, UnitRange{Int}),all=false)
    if all
        c.parameters[i,:]
    else
        c.parameters[i,c.params2s_nms]
    end
end
function parameters(c::AbstractChain;all=false)
    if all
        c.parameters
    else
        c.parameters[:,c.params2s_nms]
    end
end

# FIXME
# could do the same for moments? i.e. only return MProb.moments_subset?

moments(c::AbstractChain)                       = c.moments[i,:]
moments(c::AbstractChain, i::UnitRange{Int})    = c.moments[i,:]
moments(c::AbstractChain, i::Int)               = moments(c, i:i)
infos(c::AbstractChain)                         = c.infos
infos(c::AbstractChain, i::UnitRange{Int})      = c.infos[i,:]
infos(c::AbstractChain, i::Int)                 = infos(c, i:i)
evals(c::AbstractChain, i::UnitRange{Int})      = c.infos[i,:evals]
evals(c::AbstractChain, i::Int)                 = evals(c, i:i)
allstats(c::AbstractChain,i::UnitRange{Int})    = cbind(c.infos[i,:],c.parameters[i,:],c.moments[i,:])
allstats(c::AbstractChain,i::Int)               = cbind(infos(c, i:i),parameters(c, i:i),moments(c, i:i))

# ---------------------------  BGP GETTERS / SETTERS ----------------------------
export getEval, getLastEval, appendEval

function getEval(chain::AbstractChain, i::Int64)
    ev = Eval()
    ev.value  = chain.infos[i,:evals]
    ev.time   = chain.infos[i,:eval_time]
    ev.status = chain.infos[i,:status]

    for k in names(chain.parameters)
        if !(k in [:chain_id,:iter])
            ev.params[k] = chain.parameters[i,k]
        end
    end
    for k in names(chain.moments)
        ev.moments[k] = chain.moments[i,k]
    end

    return (ev)
end

getLastEval(chain :: AbstractChain) = getEval(chain::AbstractChain, chain.i - 1 )

function appendEval!(chain::AbstractChain, ev:: Eval, ACC::Bool, prob::Float64)
    # push!(chain.infos,[chain.i,val,ACC,status,0,prob])
    chain.infos[chain.i,:evals]  = ev.value
    chain.infos[chain.i,:prob]   = prob
    chain.infos[chain.i,:accept] = ACC
    chain.infos[chain.i,:status] = ev.status
    chain.infos[chain.i,:eval_time] = ev.time
    for im in chain.moments_nms
        chain.moments[chain.i,im] = ev.moments[im]
    end
    for ip in chain.params_nms
        chain.parameters[chain.i,ip] = ev.params[ip]
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
function parameters(MC::Array,i::Union(Integer, UnitRange{Int});allp=false)
    if !isa(MC[1],AbstractChain)
        error("must give array of AbstractChain") 
    end
    r = parameters(MC[1],i,allp) 
    if length(MC)>1
        for ix=2:length(MC)
            r = rbind(r,parameters(MC[ix],i,allp))
        end
    end
    return r
end

function parameters(MC::Array,all::Bool)
    parameters(MC,1:MC[1].i,allp=all)
end


function moments(MC::Array,i::Union(Integer, UnitRange{Int}))
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






# update the iteration count on each chain
function incrementChainIter!(MC::Array)
    for ix in 1:length(MC)
        MC[ix].i += 1
    end 
end



function saveChainToHDF5(chain::AbstractChain, ff5::HDF5File,path::ASCIIString)
    simpleDataFrameSave(chain.parameters,ff5,joinpath(path,"parameters"))
    simpleDataFrameSave(chain.infos,ff5, joinpath(path,"infos"))
    simpleDataFrameSave(chain.moments,ff5, joinpath(path,"moments"))
end

function simpleDataFrameSave(dd::DataFrame,ff5::HDF5File, path::ASCIIString)
    for nn in names(dd)
        if eltype(dd[nn]) <: Number
            write(ff5,joinpath(path,string(nn)),convert(Array{Float64,1},dd[nn])) 
        elseif eltype(dd[nn]) <: String
            write(ff5,joinpath(path,string(nn)),convert(Array{ASCIIString,1},dd[nn])) 
        end
    end
end



