
# TODO
# add plot slices


export Slice,slices,add!,get

type Slice
    res  :: Dict    # result
    p0   :: Dict    # initial param
    m0   :: Dict    # data moments

    function Slice(p::Dict,m::Dict)
        this = new()
        this.res = [ k => Dict() for k in keys(p) ]
        this.p0=deepcopy(p)
        this.m0=deepcopy(m)
        return this
    end

end

function add!(s::Slice, p::Symbol, ev::Eval)
    s.res[p][ ev.params[p] ] = Dict(:moments => ev.simMoments, :value => ev.value )
end

function get(s::Slice, p::Symbol, m::Symbol)
    x = [ k for (k,v) in s.res[p] ]

    if (m == :value)
        y =  [ v[:value] for (k,v) in s.res[p] ]
    else
        y =  [ v[:moments][m] for (k,v) in s.res[p] ]
    end

    xx = convert(Array{Float64,1},x)
    ix = sortperm(xx)

    return Dict( :x => xx[ix] , :y => convert(Array{Float64,1},y)[ix] )
end

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


#'.. py:function:: slices(m,pad)
#'
#'   computes slices for the objective function
function slices(m::MProb,npoints::Int,pad=0.1)

    t0 = time()
    ranges = m.params_to_sample

    res = Slice(m.initial_value, m.moments)
    # we want to iterate over each parameters, 
    # and compute a linerange
    for (pp,bb) in m.params_to_sample
    
        # initialize eval
        ev = Eval(m,m.initial_value)
        Lumberjack.info("slicing along $pp")

        vv = pmap( linspace(bb[:lb], bb[:ub], npoints) ) do pval
            ev2 = deepcopy(ev)
            ev2.params[pp] = pval
            ev2 = evaluateObjective(m,ev2)
            return(ev2)
        end

        for (v in vv) 
            add!( res, pp, v)
        end
    end  
    t1 = round((time()-t0)/60)
    Lumberjack.info("done after $t1 minutes")

    return res 
end

function write(s::Slice,fname::ASCIIString)

    ff5 = h5open(fname, "w") 
    # save res
    for (k,v) in s.res
        g = HDF5.g_create(ff5,string(k))
        attrs(g)["Description"] = "Values and moments along paramter $k"
        for (k2,v2) in v
            g2 = HDF5.g_create(g,string(k2))
            g2["value"] = v2[:value]
            m = HDF5.g_create(g2,"moments")
            for (km,vm) in v2[:moments]
                m[string(km)] = vm
            end
        end
    end
    close(ff5)
end

function readSlice(fname::ASCIIString)

    # get data
    ff5 = h5open(fname, "r") 

    sl = Slice(Dict("dummyp"=>0),Dict("dummym"=>2))
    pnames = names(ff5)

    res = Dict()

    for p in pnames
        g = ff5[p]
        res[symbol(p)] = Dict()
        for pp in names(g) 
            res[symbol(p)][float(pp)] = Dict()
            res[symbol(p)][float(pp)][:value] = read(g[pp]["value"])
            res[symbol(p)][float(pp)][:moments] = Dict()
            mnames = names(g[pp]["moments"])
            for mn in mnames
                md = g[pp]["moments"][mn]
                res[symbol(p)][float(pp)][:moments][symbol(mn)] = read(md)
            end
        end
    end

    sl.res = res
    return sl
end

function readSliceRemote(remote::ASCIIString) 
    a = tempname()
    run(`scp $remote $a`)
    println("saving locally to $a")
    return readSlice(a)
end
