#'.. py:function:: slices(m,pad)
#'
#'   computes slices for the objective function
function slices(m::MProb,npoints::Int,pad=0.1)

    ranges = m.params_to_sample


    res = Dict()
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

        sub_res = Dict()
        for (v in vv) 
            sub_res[:par_val] = v.params[pp]
            sub_res[:moments] = v.moments
            sub_res[:params]  = v.params
        end
        res[pp] = sub_res
    end   
   
end


