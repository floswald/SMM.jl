
# TODO
# add plot slices
type MyT
    data :: Dict
end

# @recipe function f{T<:Distribution}(::Type{T}, dist::T; func = pdf)
#     xi -> func(dist, xi)
# end






export Slice,slices,add!,get,optSlices


"""
    Slice 

A *slice* in dimension ``j`` of a function ``f \in \mathbb{R}^N`` is defined as ``f(p[1],...,p[j],...,p[N])``, where `p` is the initial parameter vector and `p[j] = linspace(lower[j],upper[j],npoints)`, where `lower[j],upper[j]` are the bounds of the parameter space in dimension ``j``.

## Fields

* `res`: `Dict` of resulting slices. For each parameter ``p_j`` there is a `Dict` with as many entries as `npoints` chosen in [`doSlices`](@ref)
* `p0`: Initial parameter vector dict
* `m0`: data moments dict

# Examples

```julia
julia> using MOpt
julia> m = MProb()
julia> p = OrderedDict(:p1=>1.1,:p2=>pi)
julia> m = OrderedDict(:m1=>rand(),:m2=>e)
julia> Slice(p,m)

```

"""
type Slice
    res  :: Union{Dict,OrderedDict}    # result
    p0   :: Union{Dict,OrderedDict}    # initial param
    m0   :: Union{Dict,OrderedDict}    # data moments

    function Slice(p::Union{Dict,OrderedDict},m::Union{Dict,OrderedDict})
        this = new()
        this.res = Dict( k => Dict() for k in keys(p) )
        this.p0=deepcopy(p)
        this.m0=deepcopy(m)
        return this
    end
end

function add!(s::Slice, p::Symbol, ev::Eval)
    s.res[p][ ev.params[p] ] = Dict(:moments => ev.simMoments, :value => ev.value )
end

function get(s::Slice, p::Symbol, m::Symbol)
    # x = [ k for (k,v) in s.res[p] ]
    x = collect(keys(s.res[p]))

    if (m == :value)
        y =  [ v[:value] for (k,v) in s.res[p] ]
    else
        y =  [ v[:moments][m] for (k,v) in s.res[p] ]
    end

    xx = convert(Array{Float64,1},x)
    ix = sortperm(xx)

    return Dict( :x => xx[ix] , :y => convert(Array{Float64,1},y)[ix] )
end

# @recipe function f(s::Slice , w::Symbol )
# # function plotit( s::Slice , w::Symbol)
#     n = length(s.res)  # num of subplots
#     z = sqrt(n)
#     rows, cols = (floor(Int,z), ceil(Int,z))
#     p = Any[]
#     for (k,v) in s.res 
#         d = get(s, k , w )
#         push!(p,plot(d[:x],d[:y],xlabel="$k",ylabel="$w",legend=false,linecolor=:black))
#     end
#     plot(p...,layout=grid(rows,cols))
# end



# MOpt.plotit(s,:value)
# function get(s::Slice, p::Symbol)
#     rr = Dict()

#     rr[:x]     = [ k for (k,v) in s.res[p] ]
#     rr[:value] = [ v[:value] for (k,v) in s.res[p] ]
#     rr[:moment]
#     for (k,v) in s.res[p][:moment]
#         rr[:moment][k] = v
#     end

#     return [ :x => convert(Array{Float64,1},x) , :y => convert(Array{Float64,1},y) ]
# end




"""
    optSlices(m::MProb,npoints::Int,parallel=false)

Computes [`Slice`](@ref)s of an [`MProb`](@ref) and keeps the best value from each slice. This implements a naive form of gradient descent in that it optimizes wrt to one direction at a time, keeping the best value. It's naive because it does a grid search in that direction. The grid size shrinks, however, at a rate `update`.
"""
function optSlices(m::MProb,npoints::Int;parallel=false,tol=1e-5,update=nothing,filename="trace.jld2")

    t0 = time()
    # res = Slice(m.initial_value, m.moments)
    # we want to iterate over each parameters, 
    # and compute a linerange

    ranges = m.params_to_sample
    cur_param = m.initial_value
    bestp = deepcopy(m.initial_value)
    dvec = deepcopy(cur_param)
    for (k,v) in dvec
        v = Inf
    end

    # output
    dout = Dict()

    # history data.frame
    df0 = DataFrame()

    delta = Inf
    iter = 0

    prog = ProgressThresh(tol, "Minimizing:")
    while delta > tol
        ProgressMeter.update!(prog, delta)
        # println("current search range:")
        # print(json(ranges,4))
        iter += 1
        println("iteration $iter")

        for (pp,bb) in ranges
            println("   working on $pp")
        
            # initialize eval
            cur_param = deepcopy(bestp)
            ev = Eval(m,cur_param)

            if parallel
                takes = @elapsed vv = pmap( linspace(bb[:lb], bb[:ub], npoints) ) do pval
                    ev2 = deepcopy(ev)
                    ev2.params[pp] = pval
                    ev2 = evaluateObjective(m,ev2)
                    return(ev2)
                end
            else
                takes = @elapsed vv = map( linspace(bb[:lb], bb[:ub], npoints) ) do pval
                    ev2 = deepcopy(ev)
                    ev2.params[pp] = pval
                    ev2 = evaluateObjective(m,ev2)
                    return(ev2)
                end
            end

            allvals = Dict()

            # find best parameter value
            minv = Inf
            bestp = deepcopy(ev.params)
            for iv in 1:length(vv)
                if (typeof(vv[iv]) <: Exception)
                        # warn("exception received. value not stored.")
                    allvals[iv] = Dict(:p => "Exception", :value => NaN)
                else
                    val = vv[iv]
                    allvals[iv] = Dict(:p => val.params, :value => val.value)
                    # println("good value? $(isfinite(val.value) && val.value < minv)")
                    if isfinite(val.value) && val.value < minv
                        minv = val.value
                        bestp = deepcopy(val.params)
                        dout[:best] = Dict(:p => val.params, :value => val.value)
                        # println("best value for $pp is $minv")
                    # else
                    #     dout[:best] = iter > 2 ? dout[:best] : Dict(:p => "Exception", :value => NaN)
                    end
                end
                # println("best value so far:")
                # print(json(dout[:best],4))
            end
            for (k,v) in allvals
                if (nrow(df0)) > 0
                    x = DataFrame(iter=iter,param=pp,val_idx=k)
                    for (ki,vi) in v[:p]
                        x[ki] = vi
                    end
                    x[:value] = v[:value]
                    append!(df0,x)
                else
                    df0[:iter] = iter
                    df0[:param] = pp
                    df0[:val_idx] = k
                    for (ki,vi) in v[:p]
                        df0[ki] = vi
                    end
                    df0[:value] = v[:value]
                end
            end
            sort!(df0,(:iter,:param,:val_idx))
            dout[:history] = df0

            if takes > 60
                JLD2.@save filename dout
            end

            dvec[pp] = cur_param[pp] - bestp[pp]
        end  # end all values in ranges

        # update search ranges
        # maintain range boundaries
        if update!=nothing
            for (k,v) in bestp
                r = (ranges[k][:ub] - ranges[k][:lb])/2
                ranges[k][:lb] = max(v - update * r,ranges[k][:lb])
                ranges[k][:ub] = min(v + update * r,ranges[k][:ub])
            end
            @debug(logger,"search ranges updated to $ranges")
        end

        # println(cur_param)
        # println(bestp)
        delta = norm(collect(values(dvec)))
    end
    t1 = round((time()-t0)/60)
    # @info(logger,"done after $t1 minutes")
    JLD2.@save filename dout
    return dout
end


"""
    doSlices(m::MProb,npoints::Int,parallel=false)

Computes [`Slice`](@ref)s of an [`MProb`](@ref)
"""
function doSlices(m::MProb,npoints::Int,parallel=false)

    t0 = time()
    res = Slice(m.initial_value, m.moments)
    # we want to iterate over each parameters, 
    # and compute a linerange
    for (pp,bb) in m.params_to_sample
    
        # initialize eval
        ev = Eval(m,m.initial_value)
        @info(logger,"slicing along $pp")

        if parallel
            vv = pmap( linspace(bb[:lb], bb[:ub], npoints) ) do pval
                ev2 = deepcopy(ev)
                ev2.params[pp] = pval
                ev2 = evaluateObjective(m,ev2)
                return(ev2)
            end
        else
            vv = map( linspace(bb[:lb], bb[:ub], npoints) ) do pval
                ev2 = deepcopy(ev)
                ev2.params[pp] = pval
                ev2 = evaluateObjective(m,ev2)
                return(ev2)
            end
        end

        for v in vv 
            if (typeof(v) <: Exception)
                warn("exception received. value not stored.")
            else
                add!( res, pp, v)
            end
        end
    end  
    t1 = round((time()-t0)/60)
    @info(logger,"done after $t1 minutes")

    return res 
end

function save(s::Slice,fname::String)
    FileIO.save(fname,"s",s)
end

function load(fname::String)
    FileIO.load(fname)
end
