
# TODO
# add plot slices
type MyT
    data :: Dict
end

m = MyT(Dict(k => 
             Dict(j => Dict(:value=>rand(),
                            :moments=>rand(3)) for j in linspace(-1,1,10))
            for k in [:a, :b])
    )

@recipe function f(::Type{MyT}, r::MyT )
function values(r::MyT )
    n = length(r.data)  # num of subplots
    z = sqrt(n)
    rows, cols = (floor(Int,z), ceil(Int,z))
    x = Dict(k => collect(keys(r.data[k])) for k in keys(r.data))

    p = Any[]
    for (k,v) in r.data 
        d = MOpt.get(r, :k , :value)
        push!(p,plot(d[:x],d[:y]))
    end
    plot(p...,layout=grid(rows,cols))
end

# @recipe function f{T<:Distribution}(::Type{T}, dist::T; func = pdf)
#     xi -> func(dist, xi)
# end






export Slice,slices,add!,get


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

"""
    doSlices(m::MProb,npoints::Int,pad=0.1)

Computes [`Slice`](@ref)s of an [`MProb`](@ref)
"""
function doSlices(m::MProb,npoints::Int,pad=0.1)

    t0 = time()
    res = Slice(m.initial_value, m.moments)
    # we want to iterate over each parameters, 
    # and compute a linerange
    for (pp,bb) in m.params_to_sample
    
        # initialize eval
        ev = Eval(m,m.initial_value)
        @info("slicing along $pp")

        vv = pmap( linspace(bb[:lb], bb[:ub], npoints) ) do pval
            ev2 = deepcopy(ev)
            ev2.params[pp] = pval
            ev2 = evaluateObjective(m,ev2)
            return(ev2)
        end

        for v in vv 
            if (typeof(v) <: Exception)
                @warn("exception received. value not stored.")
            else
                add!( res, pp, v)
            end
        end
    end  
    t1 = round((time()-t0)/60)
    @info("done after $t1 minutes")

    return res 
end

function write(s::Slice,fname::String)

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

function readSlice(fname::String)

    # get data
    ff5 = h5open(fname, "r") 

    sl = Slice(Dict("dummyp"=>0),Dict("dummym"=>2))
    pnames = names(ff5)

    res = Dict()

    for p in pnames
        g = ff5[p]
        res[Symbol(p)] = Dict()
        for pp in names(g) 
            res[Symbol(p)][float(pp)] = Dict()
            res[Symbol(p)][float(pp)][:value] = read(g[pp]["value"])
            res[Symbol(p)][float(pp)][:moments] = Dict()
            mnames = names(g[pp]["moments"])
            for mn in mnames
                md = g[pp]["moments"][mn]
                res[Symbol(p)][float(pp)][:moments][Symbol(mn)] = read(md)
            end
        end
    end

    sl.res = res
    return sl
end

function readSliceRemote(remote::String) 
    a = tempname()
    run(`scp $remote $a`)
    println("saving locally to $a")
    return readSlice(a)
end
