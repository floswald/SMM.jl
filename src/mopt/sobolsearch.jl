export sobolsearch

using Sobol

#'.. py:function:: slices(m,pad)
#'
#'   computes slices for the objective function
function sobolsearch(m::MProb,N::Int64)

    ranges = m.params_to_sample
    np     = length(keys(ranges))

    sq = SobolSeq(np)

    # get number of nodes
    nc = nprocs()
    res = Eval[]
    ev = Eval(m,m.initial_value)

    while N >=0

        SS = [ Dict(ps2s_names(m),next(sq)) for i in 1:(2*nc)]
        vv = pmap(SS) do ushock
            ev2 = deepcopy(ev)

            for (param,bound) in ranges
                ev2.params[param] = bound[:lb] + ushock[param]*(bound[:ub] - bound[:lb]) 
            end

            ev2 = evaluateObjective(m,ev2)
            return(ev2)
        end

        N-=length(SS)
        append!(res,vv)
    end

    return(res)   
end


