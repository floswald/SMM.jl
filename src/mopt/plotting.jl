


# @shorthands slices

# @recipe function f(::Type{Val{:slices}}, plt::Plot; density = false)
    

#     # set up the subplots
#     legend --> false
#     link := :both
#     grid --> false
#     layout --> @layout [
#         tophist           _
#         hist2d{0.9w,0.9h} righthist
#     ]

#     # main histogram2d
#     @series begin
#         seriestype := :histogram2d
#         right_margin --> 0mm
#         top_margin --> 0mm
#         subplot := 2
#     end

#     # these are common to both marginal histograms
#     ticks := nothing
#     guide := ""
#     foreground_color_border := nothing
#     fillcolor --> Plots.fg_color(d)
#     linecolor --> Plots.fg_color(d)

#     if density
#         trim := true
#         seriestype := :density
#     else
#         seriestype := :histogram
#     end

#     # upper histogram
#     @series begin
#         subplot := 1
#         bottom_margin --> 0mm
#         y := x
#     end

#     # right histogram
#     @series begin
#         orientation := :h
#         subplot := 3
#         left_margin --> 0mm
#         y := y
#     end
# end

# # now you can plot like:
# marginalhist(rand(1000), rand(1000))

function best_grid(n)
    rows = floor(Int,sqrt(n))
    cols = ceil(Int,n/rows)
    rows,cols
end

# @recipe function f{T<:Distribution}(::Type{T}, dist::T; func = pdf)
#     xi -> func(dist, xi)
# end



@recipe function f(ma::MAlgoBGP,chain::Int)
    # produce a subplot for each parameter
    n = length(ma.m.params_to_sample)  # num of subplots
    rows,cols = best_grid(n)
    mat = zeros(Int,rows,cols)
    empty = mod(length(mat),n)

    # build a grid
    g = grid(rows,cols, heights=ones(rows)/rows,widths=ones(cols)/cols)
    if empty > 0
        idx = 1
        for i=1:rows,j=1:cols
            if idx > length(mat) - empty
                g[i,j].attr[:blank] = true
            end
            idx += 1
        end
    end
    spidx = 1
    indices = Dict()
    for (k,v) in ma.m.params_to_sample
        indices[k] = spidx
        spidx += 1
    end
    layout := g    
    # foreground_color_border := nothing
    margin    --> 1mm
    titlefont --> font(11)
    fillcolor --> :orange
    fillrange --> 0
    fillalpha --> 0.5
    grid      --> false 
    link      --> :none
    xticks     := true
    xguidefont  := font(10)
    legend := false

    dd = MOpt.cur_param(ma)[chain]
    for (k,v) in ma.m.params_to_sample
        dist = Normal(dd[:mu][k],diag(dd[:sigma])[indices[k]])
        println(dist)
        x = linspace(v[:lb],v[:ub],100)
        @series begin
            seriestype := :path
            subplot := indices[k]
            x,pdf(dist,x)
        end
        @series begin
            subplot := indices[k]
            seriestype := :line
            linetype := :vline
            linecolor := :red
            xguide     := "$k"
            [dd[:mu][k]]
        end
    end

end


@recipe function f(s::Slice , w::Symbol )
    n = length(s.res)  # num of subplots
    rows,cols = best_grid(n)
    mat = zeros(Int,rows,cols)
    empty = mod(length(mat),n)


    # build a grid
    g = grid(rows,cols, heights=ones(rows)/rows,widths=ones(cols)/cols)
    if empty > 0
        idx = 1
        for i=1:rows,j=1:cols
            if idx > length(mat) - empty
                g[i,j].attr[:blank] = true
            end
            idx += 1
        end
    end

    spidx = 1
    indices = Dict()
    for (k,v) in s.res 
        indices[k] = spidx
        spidx += 1
    end
    layout := g

    # labels
    # labs = pop!(d, :label, [""])  # d is a special KW Dict

    # some defaults
    legend    := false
    # foreground_color_border := nothing
    margin    --> 1mm
    titlefont --> font(11)
    fillcolor --> Plots.fg_color(d)
    linecolor --> Plots.fg_color(d)
    grid      --> true
    link      --> :x
    xticks     := true
    xguidefont  := font(10)

    # build series
    for (k,v) in s.res 
        da = get(s,k,w)
        @series begin
            seriestype := :line
            subplot    := indices[k]
            xguide     := "$k"
            # xticks     := pos[k][1]==rows ? true : false
            # xguide     := pos[k][1]==rows ? "$k" : nothing
            yguide     := mod(indices[k],cols)==1 ? "$w" : ""
            da[:x], da[:y]
        end
        if haskey(d,:markersize)
            @series begin
                seriestype := :scatter
                subplot    := indices[k]
                xguide     := "$k"
                # xticks     := pos[k][1]==rows ? true : false
                # xguide     := pos[k][1]==rows ? "$k" : nothing
                yguide     := mod(indices[k],cols)==1 ? "$w" : ""
                da[:x], da[:y]
            end
        end
    end
end

function testslice()

s = Slice(Dict(:p1=>1, :p2=>0.0),Dict(:m1=>rand(),:m2=>rand()))
s.res = Dict(k => Dict(j => Dict(:value=>rand(),:moments=>Dict(:m1=>rand(),:m2=>rand())) for j in linspace(-1,1,10)) for k in [Symbol("p$i") for i in 1:17])
MOpt.plot(s,:value)
# MOpt.plot(s,:value,markersize=10,markercolor=:red,link=:both)
end

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


function plotSlices(m::MProb,val_df::DataFrame,mom_df::DataFrame,facet::Array{String,1})

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
plotSlices(m::MProb,val_df::DataFrame,mom_df::DataFrame,facet::String) = plotSlice(m,val_df,mom_df,[facet])

