

mutable struct MyT
    data :: Dict
end

# @recipe function f{T<:Distribution}(::Type{T}, dist::T; func = pdf)
#     xi -> func(dist, xi)
# end






export Slice,slices,add!,get


"""
    Slice 

A *slice* in dimension ``j`` of a function ``f ∈ \\mathbb{R}^N`` is defined as ``f(p[1],...,p[j],...,p[N])``, where `p` is the initial parameter vector and `p[j] = linspace(lower[j],upper[j],npoints)`, where `lower[j],upper[j]` are the bounds of the parameter space in dimension ``j``.

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
mutable struct Slice
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
    doSlices(m::MProb,npoints::Int,parallel=false,pad=0.1)

Computes [`Slice`](@ref)s of an [`MProb`](@ref)
"""
function doSlices(m::MProb,npoints::Int,parallel=false,pad=0.1)

    t0 = time()
    res = Slice(m.initial_value, m.moments)
    # we want to iterate over each parameters, 
    # and compute a linerange
    for (pp,bb) in m.params_to_sample
    
        # initialize eval
        ev = Eval(m,m.initial_value)
        @info "slicing along $pp"

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
                @warn "exception received. value not stored."
            else
                add!( res, pp, v)
            end
        end
    end  
    t1 = round((time()-t0)/60)
    @info "done after $t1 minutes"

    return res 
end

function save(s::Slice,fname::String)
    FileIO.save(fname,"s",s)
end

function load(fname::String)
    FileIO.load(fname)
end
