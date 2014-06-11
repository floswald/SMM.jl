
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
function collectFields(dict::Dict, I::UnitRange{Int}, df::Bool=false)
    if length(I) == 0
        println("no evaluations to show")
    end
    if df
        cols = [dict[k][I] for k in keys(dict)]
        cnames = Array(Symbol,length(dict))
        pkeys = collect(keys(dict))
    for i in 1:length(pkeys)
        cnames[i] = symbol(pkeys[i])
    end
        return DataFrame(cols, cnames)
    else ## ==== return as collection
        return({ k => v[I] for (k,v) in dict })
    end
end

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





## Multiple Default Chains
## =======================

# Stores multilpe chains
type MChain
  n :: Int # number of chains
  chains :: Array

  # function MChain(n,ChType::AbstractChain,MProb,L)
  #   chains = [ ChType(MProb,L) for i in 1:n ]
  #   return new(n,chains)
  # end
  function MChain(n,MProb,L)
    chains = [ Chain(MProb,L) for i in 1:n ]
    return new(n,chains)
  end
end

# can also return a range of the dataframe
function getindex(mc::MChain, i::Int)
    return mc.chains[i]
end

# methods for MChain
function appendEval!(MC::MChain, which::Int, vals::Dict, acc::Bool,status::Int)
    appendEval!(MC.chains[which],vals,acc,status)
end

# update the iteration count on each chain
function updateIter!(MC::MChain)
    for ix in 1:MC.n
        MC.chains[ix].i += 1
    end 
end

