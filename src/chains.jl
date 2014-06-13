
# the goal here is to have a convenient way of 
# storing the different evaluations.
# That includes the realizations of the objective
# function, the moments, and the value of the parameters
# that went in the evaluation.
# It might also include additional information
#
# so usually we'll have to refer to a chain number, an iteration, and an object
# which might be either a moment, a parameter, or an output.
#
# Because the chains might grow over time and because storing a in a DataFrame 
# is not always easy, we are going to use vectors for each.


# defining an abstract chain type in case 
# some algorithm need additional information
# on top of (eval, moment, parameter)
abstract AbstractChain

# methods for AbstractChain
# -------------------------

# taking a dictionary of vectors, returns
# the values as a dataframe or as a dictionary
# using Debug
# @debug function collectFields(dict::Dict, I::UnitRange{Int}, df::Bool=false)
function collectFields(dict::Dict, I::UnitRange{Int}, df::Bool=false)
    # if length(I) == 0
    #     println("no evaluations to show")
    # end
    n = length(dict)
    dk = collect(keys(dict))
    if df
        cols = Any[dict[k][I] for k in dk]  # notice: Any is crucial here to get type-stable var
        # cols = Array(Any,n)
        # for i in 1:n
        #     cols[i] = dict[dk[i]]
        # end
        cnames = Array(Symbol,length(dict))
        for i in 1:n
            cnames[i] = symbol(dk[i])
        end
        return DataFrame(cols, cnames)
    else ## ==== return as collection
        return({ k => v[I] for (k,v) in dict })
    end
end

# taking a dataframe row
# fills in the values into keys of a dict at row I of 
# arrays in dict
function fillinFields!(dict::Dict,df::DataFrame,I::Int)

    if nrow(df)!=1
        error("can fill in only a single dataframe row")
    end
    dk = collect(keys(dict))
    for ik in dk
        dict[ik][I] = df[symbol(ik)][1]
    end

end

# same but for dict with only on entry per key
function fillinFields!(dict::Dict,df::DataFrame)

    if nrow(df)!=1
        error("can fill in only a single dataframe row")
    end
    dk = collect(keys(dict))
    for ik in dk
        dict[ik] = df[symbol(ik)][1]
    end

end

# dataframe to dict function
function df2dict(df::DataFrame)
  nm = names(df)
  snm = map(x->string(x),nm)
  out ={i => df[symbol(i)] for i in snm}
  return out
end


# methods for a single chain
# ==========================

collectFields(dict::Dict,df::Bool=false)                 = collectFields(dict, 1:length(dict),df)
collectFields(dict::Dict,i::Int, df::Bool=false)         = collectFields(dict, i:i, df)
parameters(c::AbstractChain, df::Bool=false)                     = collectFields(c.parameters, 1:c.i, df)
parameters(c::AbstractChain, I::UnitRange{Int}, df::Bool=false)  = collectFields(c.parameters, I, df)
parameters(c::AbstractChain, i::Int, df::Bool=false)             = collectFields(c.parameters, i, df)
moments(c::AbstractChain, df::Bool=false)                        = collectFields(c.moments, 1:c.i, df)
moments(c::AbstractChain, I::UnitRange{Int}, df::Bool=false)     = collectFields(c.moments, I, df)
moments(c::AbstractChain, i::Int, df::Bool=false)                = collectFields(c.moments, i, df)
infos(c::AbstractChain, df::Bool=false)                          = collectFields(c.infos, 1:c.i, df)
infos(c::AbstractChain, I::UnitRange{Int}, df::Bool=false)       = collectFields(c.infos, I, df)
infos(c::AbstractChain, i::Int, df::Bool=false)                  = collectFields(c.infos, i, df)
alls(c::AbstractChain, df::Bool=false)                           = collectFields(merge(c.infos,c.parameters,c.moments), 1:c.i, df)
alls(c::AbstractChain, I::UnitRange{Int}, df::Bool=false)        = collectFields(merge(c.infos,c.parameters,c.moments), I, df)
alls(c::AbstractChain, i::Int, df::Bool=false)                   = collectFields(merge(c.infos,c.parameters,c.moments), i, df)
getindex(c::AbstractChain, i::UnitRange{Int}) = alls(c,i,true)
getindex(c::AbstractChain, i::Int) = alls(c,i,true)
evals(c::AbstractChain, df::Bool=false)                          = collectFields({"evals" => c.infos["evals"]}, 1:c.i, df)["evals"]
evals(c::AbstractChain, I::UnitRange{Int}, df::Bool=false)       = collectFields({"evals" => c.infos["evals"]}, I, df)["evals"]
evals(c::AbstractChain, i::Int, df::Bool=false)                  = collectFields({"evals" => c.infos["evals"]}, i, df)["evals"]
status(c::AbstractChain, df::Bool=false)                          = collectFields({"status" => c.infos["status"]}, 1:c.i, df)["status"]
status(c::AbstractChain, I::UnitRange{Int}, df::Bool=false)       = collectFields({"status" => c.infos["status"]}, I, df)["status"]
status(c::AbstractChain, i::Int, df::Bool=false)                  = collectFields({"status" => c.infos["status"]}, i, df)["status"]
accept(c::AbstractChain, df::Bool=false)                          = collectFields({"accept" => c.infos["accept"]}, 1:c.i, df)["accept"]
accept(c::AbstractChain, I::UnitRange{Int}, df::Bool=false)       = collectFields({"accept" => c.infos["accept"]}, I, df)["accept"]
accept(c::AbstractChain, i::Int, df::Bool=false)                  = collectFields({"accept" => c.infos["accept"]}, i, df)["accept"]

# appends values from objective function
# at CURRENT iteration
function appendEval!(chain::AbstractChain, vals::Dict, ACC::Bool, status::Int)
  chain.infos["evals"][chain.i]  = vals["value"]
  chain.infos["accept"][chain.i] = ACC
  chain.infos["status"][chain.i] = status
  for (k,v) in vals["moments"]
    chain.moments[k][chain.i] = v
  end
  for (k,v) in vals["params"]
    chain.parameters[k][chain.i] = v
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
function Allparameters(MC)
    # TODO how to check that MC is an array of AbstractChains?
    # if super(typeof(MC[1]))!=AbstractChain 
    #     error("must give array of AbstractChain") 
    # end
    r0 = parameters(MC[1],true) # collects all params up to current iteration i
    r = cbind(DataFrame(id=[1 for i=1:nrow(r0)],iter=1:nrow(r0)),r0)
    if length(MC)>1
        for ix=2:length(MC)
            r = rbind(r,cbind(DataFrame(id=[ix for i=1:nrow(r0)],iter=1:nrow(r0)),parameters(MC[ix],true)))
        end
    end
    return r
end


# the default chain type
# we create a dictionary with arrays
# for each parameters
type Chain <: AbstractChain
  i::Int              # current index
  infos      ::Dict   # dictionary of arrays(L,1) with eval, ACC and others
  parameters ::Dict   # dictionary of arrays(L,1), 1 for each parameter
  moments    ::Dict   # dictionary of DataArrays(L,1), 1 for each moment

  function Chain(MProb,L)
    infos      = { "evals" => @data([0.0 for i = 1:L]) , "accept" => @data([false for i = 1:L]), "status" => [0 for i = 1:L]}
    parameters = { x => zeros(L) for x in ps_names(MProb) }
    moments    = { x => @data([0.0 for i = 1:L]) for x in ms_names(MProb) }
    return new(0,infos,parameters,moments)
  end
end






# trying to understand parametric types:
# abstract MyAbstract

# type MySubType <: MyAbstract
#     field1 :: ASCIIString
#     function MySubType(x)
#         new(x)
#     end
# end


# type Myt{ T<:MyAbstract}
#     chain :: Array{T,1}
#     function Myt(n,z::T) 
#         x = [z for i=1:n]
#         new(x)
#     end
# end
# Myt{T}(n::Integer,z::T) = Myt{T}(n,z)

# function newType(whichtype,x)
#     if isa(whichtype,MySubType)
#         return MySubType(x)
#     else
#         println("no other types")
#     end
# end


# this doesnt work
# but really we can just use Array{Chain,1}

## Multiple Default Chains
## =======================

# Stores multilpe chains
# type MChain{ T <: AbstractChain}
#   n :: Int # number of chains
#   chains :: Array{T,1}

#   function MChain(  n,ch::T,MProb,L)
#     chains = [ ch for i in 1:n ]
#     return new(n,chains)
#   end
#   # function MChain(n,MProb,L)
#   #   chains = [ Chain(MProb,L) for i in 1:n ]
#   #   return new(n,chains)
#   # end
# end
# MChain{T}(n::Integer,ch::T,MProb::MProb,L::Integer) = MChain{T}(n,ch,MProb,L)


# function newChain(whichChain,m::MProb,L::Integer)
#     if isa(whichChain,Chain)
#         return Chain(m,L)
#     else
#         println("no other chain types")
#     end
# end

# can also return a range of the dataframe
function getindex(mc::Array{AbstractChain}, i::Int)
    return mc[i]
end

# methods for MChain
# function appendEval!(MC::Array{AbstractChain,1}, which::Int, vals::Dict, acc::Bool,status::Int)
#     appendEval!(MC[which],vals,acc,status)
# end

# update the iteration count on each chain
function updateIter!(MC::Array)
    for ix in 1:length(MC)
        MC[ix].i += 1
    end 
end




# using DataFrames

# function collectFields(dict::Dict)
#     di_keys = collect(keys(dict))
#     n = length(dict)
#     cols = {dict[k] for k in di_keys}
#     # cols = Array(Any,n)
#     # for i in 1:n
#     #     cols[i] = dict[di_keys[i]]
#     # end
#     cnames = Array(Symbol,length(dict))
#     for i in 1:length(di_keys)
#         cnames[i] = symbol(di_keys[i])
#     end
#     # return DataFrame(cols, DataFrames.Index(cnames))
#     return DataFrame(cols, cnames)
# end

# di = {ASCIIString,Array{Real,1})["a"=>[1,3],"b"=>[0.0,1.0]]
# collectFields(di)

# di2 = ["a"=>[1,3],"b"=>[0,1]]
# collectFields(di2)

