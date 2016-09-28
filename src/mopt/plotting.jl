

function plotExample(s::Slice)

    subplot(231)
    r = MOpt.get(sl, :m1 , :m1); PyPlot.plot(r[:x],r[:y],".")
    subplot(232)
    r = MOpt.get(sl, :m1 , :m2); PyPlot.plot(r[:x],r[:y],".")
    subplot(233)
    r = MOpt.get(sl, :m1 , :value); PyPlot.plot(r[:x],r[:y],".")
    subplot(234)
    r = MOpt.get(sl, :m2 , :m1); PyPlot.plot(r[:x],r[:y],".")
    subplot(235)
    r = MOpt.get(sl, :m2 , :m2); PyPlot.plot(r[:x],r[:y],".")
    subplot(236)
    r = MOpt.get(sl, :m2 , :value); PyPlot.plot(r[:x],r[:y],".")


end


function plotSlices(m::MProb,val_df::DataFrame,mom_df::DataFrame,facet::Array{ASCIIString,1})

    if "moments" in facet

        # each param has it's own figure
        # within each figure, there are nmoms subplots
        nmoms= length(m.moments_subset)
        nrows = floor(sqrt(nmoms))
        ncols = ceil(nmoms/nrows)

        # plots of moments vs parameter values
        for m_subdf in groupby(mom_df,:p_name)
            figure()
            pid = 0
            for p_subdf in groupby(m_subdf,:m_name)
                pid += 1
                subplot(nrows,ncols,pid)
                plot(p_subdf[:p_val],p_subdf[:m_val])
                axhline(y=m.moments[m.moments[:moment].==p_subdf[1,:m_name],:data_value],color="r",linestyle="--")
                axvline(x=m.initial_value[p_subdf[1,:p_name]],linestyle="-",color="grey")
                title(p_subdf[1,:m_name])
            end
            suptitle("Parameter: $(m_subdf[1,:p_name]) vs Moments")
        end

    elseif "objective" in facet

        # plots of objective function vs parameter value
        figure()
        pid = 0
        for p_subdf in groupby(val_df,:p_name)
            pid += 1
            subplot(nrows,ncols,pid)
            plot(p_subdf[:p_val],p_subdf[:f_val])
            axvline(x=m.initial_value[p_subdf[1,:p_name]],linestyle="-",color="grey")
            title(p_subdf[1,:p_name])
        end
        suptitle("Objective Function vs Parameters")

    elseif "params" in facet
        # each moment has it's own figure
        # within each figure, there are npars subplots
        npars = length(m.p2sample_sym)
        nrows = floor(sqrt(npars))
        ncols = ceil(npars/nrows)

        # plots of moments vs parameter values
        for m_subdf in groupby(mom_df,:m_name)
            figure()
            pid = 0
            for p_subdf in groupby(m_subdf,:p_name)
                pid += 1
                subplot(nrows,ncols,pid)
                plot(p_subdf[:p_val],p_subdf[:m_val])
                axhline(y=m.moments[m.moments[:moment].==p_subdf[1,:m_name],:data_value],color="r",linestyle="--")
                axvline(x=m.initial_value[p_subdf[1,:p_name]],linestyle="-",color="grey")
                title(p_subdf[1,:p_name])
            end
            suptitle("Moment: $(m_subdf[1,:m_name]) vs Parameters")
        end

        # plots of objective function vs parameter value
        figure()
        pid = 0
        for p_subdf in groupby(val_df,:p_name)
            pid += 1
            subplot(nrows,ncols,pid)
            plot(p_subdf[:p_val],p_subdf[:f_val])
            axvline(x=m.initial_value[p_subdf[1,:p_name]],linestyle="-",color="grey")
            title(p_subdf[1,:p_name])
        end
        suptitle("Objective Function vs Parameters")

    else
        ArgumentError("facet needs to be either params or moments")
    end
end

plotSlices(m::MProb,val_df::DataFrame,mom_df::DataFrame) = plotSlice(m,val_df,mom_df,"moments")
plotSlices(m::MProb,val_df::DataFrame,mom_df::DataFrame,facet::ASCIIString) = plotSlice(m,val_df,mom_df,[facet])

