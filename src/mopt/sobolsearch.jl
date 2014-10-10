export sobolsearch

using Sobol

#'.. py:function:: slices(m,pad)
#'
#'   computes slices for the objective function
function sobolsearch(m::MProb,Ntot::Int64)

    ranges = m.params_to_sample
    np     = length(keys(ranges))

    sq = SobolSeq(np)

    # get number of nodes
    nc = nprocs()
    res = Eval[]
    ev = Eval(m,m.initial_value)

    ev.value = Inf
    best_ev = ev

    N = Ntot
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

        # showing some info
        for e in vv
            if (e.status>0) & (e.value < best.value)
                best=e
            end
        end

        info(" iter left:$(Ntot - N) bestval:$(best.value) bestpar:$(best.params)")

    end

    return(res)   
end


