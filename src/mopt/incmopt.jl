

# transpose a 2-column dataframe to a one-row dataframe, so that 
# col x become new column names
# function transpose(x::DataFrame,newNames::Int)
#     if ncol(x) != 2
#         throw(ArgumentError("x must have 2 columns"))
#     end
#     if !in(newNames,[1,2])
#         throw(ArgumentError("newNames indexes col with new names: either 1 or 2"))
#     end
#     newrows = setdiff([1,2],newNames)[1]
#     # make new column names out of col index newNames
#     z = DataFrame(Float64,1,nrow(x))
#     names!(z,Symbol[y for y in x[:,newNames]])
#     for i in 1:nrow(x)
#         z[1,i] = x[i,newrows]
#     end
#     return z
# end


# taking a dictionary of vectors, returns
# the values as a dataframe or as a dictionary
# using Debug
# @debug function collectFields(dict::Dict, I::UnitRange{Int}, df::Bool=false)
function collectFields(dict::Dict, I::UnitRange{Int}, df::Bool=false)
    # if length(I) == 0
    #     println("no evaluations to show")
    # end
    n = length(dict)
    dk = sort(collect(keys(dict)))
    if df
        cols = Any[dict[k][I] for k in dk]  # notice: Any is crucial here to get type-stable var
        # cols = Array(Any,n)
        # for i in 1:n
        #     cols[i] = dict[dk[i]]
        # end
        cnames = Symbol[x for x in dk]
        return DataFrame(cols, cnames)
    else ## ==== return as collection
        return(Dict( k => v[I] for (k,v) in dict ))
    end
end

# taking a dataframe row
# fills in the values into keys of a dict at row I of 
# arrays in dict
function fillinFields!(dict::Dict,df::DataFrame,I::Int)

    if nrow(df)!=1
        ArgumentError("can fill in only a single dataframe row")
    end
    dk = collect(keys(dict))
    for ik in dk
        dict[ik][I] = df[symbol(ik)][1]
    end

end

# same but for dict with only on entry per key
# looks for keys(dict) in rownames of dataframe
function fillinFields!(dict::Dict,df::DataFrame)

    if nrow(df)!=1
        ArgumentError("can fill in only a single dataframe row")
    end
    dk = names(df)
    for ik in dk
        dict[string(ik)] = df[ik][1]
    end

end

# dataframe to dict function
function df2dict(df::DataFrame)
  nm = names(df)
  snm = map(x->string(x),nm)
  out = Dict(i => df[symbol(i)] for i in snm)
  return out
end

"""
    initsd(upper::Float64,init::Float64;prob::Float64=0.975)

finds standard deviation `sigma` of N(init,sigma) such that `prob` of the resulting pdf lies within `[init-upper,init+upper]`. This assumes a **symmetric** box around `init`.
"""
function initsd(upper::Float64,init::Float64;prob::Float64=0.975)
    @assert (prob < 1) && (prob > 0)
    @assert upper > init
    (upper - init) / quantile(Normal(),prob)
end

"""
    initvar(bound::Float64,init::Float64,prob::Float64;low_var::Float64=0.01,hi_var::Float64=100.0)

finds standard deviation `sigma` of mean zero univariate normal such that `cdf(Normal(init,sigma),bound) = prob`. 
Use this function to specify the variance when sampling from a symmetrically bounded normal distribution when you want 
to have a mass of `prob` at the bounds.
"""
function initvar(bound::Float64,init::Float64,prob::Float64;low_var::Float64=eps(),hi_var::Float64=100.0)
    @assert (prob < 1) && (prob > 0)
    fzero(x->distvar(x,init,bound,prob),low_var,hi_var)
end

"""
   distvar(sig::Float64,p::Float64)  

helper function for `initvar`
"""
function distvar(sig::Float64,mu::Float64,b::Float64,p::Float64)
    N = Normal(mu,sig)
    cdf(N,b) - p
end



# function checkbounds!(df::DataFrame,di::Dict)
# 	if nrow(df) > 1
# 		error("can only process a single row")
# 	end
# 	dfbounds = collectFields(di,1:length(di),true)
# 	for c in names(df)
# 		if df[1,c] > dfbounds[2,c]
# 			df[1,c] = dfbounds[2,c]
# 		elseif df[1,c] < dfbounds[1,c]
# 			df[1,c] = dfbounds[1,c]
# 		end
# 	end
# end

function fitMirror(x,lb,ub)
    if (x > ub)
        x2 = ub - mod( x - ub, ub - lb)
    elseif (x < lb)
        x2 = lb + mod( lb - x, ub - lb)
    else
        x2 = x
    end
    return x2 
end

fitMirror(x::Float64,d::Dict) = fitMirror(x,d[:lb],d[:ub])

function fitMirror!(x::DataFrame,b::DataFrame)
    for i in 1:length(x)
        x[i] = fitMirror(convert(Float64,x[i][1]),convert(Float64,b[i,:lb][1]),convert(Float64,b[i,:ub][1]))
    end
end

function fitMirror!(x::DataFrame,b::Dict)
    for col in names(x)
        x[1,col] = fitMirror( 
                convert(Float64,x[1,col]),
                convert(Float64,b[col][:lb]),
                convert(Float64,b[col][:ub]))
    end
end


function findInterval(x::T, vec::Array{T}) where T <: Number

    out = zeros(Int,length(x))
    sort!(vec)

    for j in 1:length(x)
        if x[j] < vec[1]
            out[1] = 0
        elseif x[j] > vec[end]
            out[end] = 0
        else
            out[j] = searchsortedfirst(vec,x[j])-1 
        end
    end
    return out
end

            

