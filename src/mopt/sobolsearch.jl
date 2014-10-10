export sobolsearch

using Sobol

#'.. py:function:: sobolsearch(m,pad)
#'
#'   This is an helping function that just evaluates the
#'   the objective function over a Sobol sequence in the 
#'   multidimensional space of parameters. This can be useful
#'   as a first pass.
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
            if (e.status>0) & (e.value < best_ev.value)
                best_ev=e
            end
        end

        info( string(@sprintf( "iter:%3i bestval:%4.2f " , Ntot - N, best_ev.value), "bestpar:$(best_ev.params)" ))

    end

    return(res)   
end


