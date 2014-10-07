
export Slice,slices,add!,get

type Slice
    res  :: Dict
    p0   :: Dict
    m0   :: Dict

    function Slice(p::Dict,m::Dict)
        this = new()
        this.res = [ k => Dict() for k in keys(p) ]
        this.p0=deepcopy(p)
        this.m0=deepcopy(m)
        return this
    end
end

function add!(s::Slice, p::Symbol, ev::Eval)
    s.res[p][ ev.params[p] ] = [:moments => ev.moments, :value => ev.value ]
end

function get(s::Slice, p::Symbol, m::Symbol)
    x = [ k for (k,v) in s.res[p] ]

    if (m == :value)
        y =  [ v[:value] for (k,v) in s.res[p] ]
    else
        y =  [ v[:moments][m] for (k,v) in s.res[p] ]
    end

    return [ :x => convert(Array{Float64,1},x) , :y => convert(Array{Float64,1},y) ]
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

    ranges = m.params_to_sample

    res = Slice(m.initial_value, m.moments)
    # we want to iterate over each parameters, 
    # and compute a linerange
    for (pp,bb) in m.params_to_sample
    
        # initialize eval
        ev = Eval(m,m.initial_value)

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

    return res 
   
end


