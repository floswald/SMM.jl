#'.. py:function:: slices(m,pad)
#'
#'   computes slices for the objective function
function slices(m::MProb,npoints::Int,pad=0.1)

    # make a dict of grids for each param
    # loop over params!
    pranges = Dict{Symbol,Array{Float64,1}}()
    for bd in m.params_to_sample
        lb = irow[:lb][1]
        ub = irow[:lb][1]
        pranges[irow[:param]] = linspace(irow[:lb][1], irow[:ub][1], npoints)  
    end
    # return pranges
    val_df = DataFrame()
    mom_df = DataFrame()
    for (k,v) in pranges
        println("currently computing slices over $k")
        dtmp = computeSlice(m,k,v)
        val_df = rbind(val_df,dtmp[1])
        mom_df = rbind(mom_df,dtmp[2])
    end
    return (val_df,mom_df)
   
end


#'.. py:function:: computeSlice(m,par,prange)
#'
#'   computes slices for the objective function
function computeSlice(m::MProb,par::Symbol,prange::Array{Float64,1})

    npar = length(prange)
    nmom = length(m.moments_subset)

    # create a list of parameter sets that cover prange
    pp = map(1:npar) do i
        p=deepcopy(m.initial_value);
        p[par]=prange[i] 
    end
        
    # evaluate the parameters in parallel
    v = pmap(x -> evaluateObjective(m,x), pp)

    # collect the different values
    mom_df = DataFrame(p_name = [par for i=1:(npar*nmom)], m_name=ASCIIString[ i for j=1:npar, i in m.moments_subset][:], p_val = repeat(prange,inner=[1],outer=[nmom]), m_val = zeros(npar*nmom))    
    val_df = DataFrame(p_name = [par for i=1:(npar)], p_val = prange, f_val = zeros(npar), status = zeros(npar))

    for ip in 1:npar
        # fill in function values
        val_df[ip, :f_val ] = v[ip]["value"]
        val_df[ip, :status] = v[ip]["status"]

        for im in 1:nmom
            # fill in moments values
            mom_df[(mom_df[:p_val].==prange[ip]) & (mom_df[:m_name].==m.moments_subset[im]), :m_val ] = v[ip]["moments"][1,symbol(m.moments_subset[im])]
        end
    end

    return (val_df,mom_df)

end
