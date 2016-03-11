export sobolsearch,sobolWeightedSearch

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

# this returns a value associated with 
function rateNewPoint(x,X,V,ni)

    if (ni<=0) 
        return 1.0
    end

    # find the nearest point
    best_val = Inf
    best_idx = 0
    for i in 1:ni

      if (   (best_val > sum((x - vec(X[i,:])).^2)) & isfinite(V[i,1]) )
        best_val = sum((x - vec(X[i,:])).^2)
        best_idx = i
      end
    end

    # returns distance times level
    return  (best_val) / (V[best_idx,1])
end


function sobolWeightedSearch(m::MProb,Ntot::Int64)
    # I want to extend Sobol to decide strategically
    # on whether to evaluate a given point based on a 
    # a tradeoff between discrepancy and levels

    ranges = m.params_to_sample
    np     = length(keys(ranges))
    sq     = SobolSeq(np)
    nc     = nprocs()
    CK     = 2*nc # chcunk size

    # get number of nodes
    res = Eval[]
    ev = Eval(m,m.initial_value)

    ev.value = Inf
    best_ev = ev

    N = Ntot
    X = zeros(N,np)
    V = zeros(N,3)
    cur_dens = 0

    while N >= CK
        ni = Ntot - N +1

        # picking new points to evaluate - this part is key
        # for each new point we want to compute a local density
        # and a local level. 
        SS = Any[]
        i  = 1
        while (i < CK)
            candidate = next(sq)
            r = rateNewPoint(candidate,X,V,ni-1)

            # compute accept reject based on average wdensity
            if ( true) #  rand() < (r/cur_dens) )
                u = next(sq)
                push!(SS, Dict(ps2s_names(m),u) )
                # update density
                cur_dens = 0.9 * cur_dens + 0.1*r
                V[ (ni+i-1) , 2] = r
                V[ (ni+i-1) , 3] = cur_dens
                X[ (ni+i-1) , :] = u
                i = i+1
            end
        end

        # map suniform shock into param space
        # and evaluate in parallel
        vv = pmap(SS) do ushock
            ev2 = deepcopy(ev)

            for (param,bound) in ranges
                ev2.params[param] = bound[:lb] + ushock[param]*(bound[:ub] - bound[:lb]) 
            end

            ev2 = evaluateObjective(m,ev2)
            return(ev2)
        end

        # update vector of values
        i = 1
        for ev in vv
            V[ (ni+i-1) ,1] = ev.value
            i = i+1
        end

        N-=length(SS)
        append!(res,vv)

        # showing some info
        for e in vv
            if (e.status>0) & (e.value < best_ev.value)
                best_ev=e
            end
        end

        println( string(@sprintf( "iter:%3i bestval:%4.2f dens:%f" , Ntot - N, best_ev.value, cur_dens ), "bestpar:$(best_ev.params)" ))

    end

    return Dict( :evals => res, :X => X, :V=>V)
end


