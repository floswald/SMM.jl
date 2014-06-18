
# define a Test objective function
function Testobj(x::Dict,mom::Dict,whichmom::Array{ASCIIString,1})

	t0 = time()
	mm = DataFrame(name=collect(keys(mom)),data=map(x -> x[1],values(mom)),model=[x["a"] + x["b"] + i for i=1:length(mom)])

	# get model moments as a dict
	mdict = {i => mm[findin(mm[:name],[i]),:model][1] for i in collect(keys(mom)) }

	# subset to required moments only
	mm = mm[findin(mm[:name],whichmom),:]

	# compute distance
	v = sum((mm[:data] - mm[:model]).^2)

	# status
	status = 1

	# time out
	t0 = time() - t0

	# return a dict
	ret = ["value" => v, "params" => x, "time" => t0, "status" => status, "moments" => mdict]
	return ret

end

 # ret = ["value" => 1.1, "params" => ["a"=>1.1,"b"=>12.1], "time" => 0, "status" => 1, "moments" => ["alpha"=>1.1,"beta"=>12.1,"gamma"=>12.1] ]


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
        return({ k => v[I] for (k,v) in dict })
    end
end

# taking a dataframe row
# fills in the values into keys of a dict at row I of 
# arrays in dict
function fillinFields!(dict::Dict,df::DataFrame,I::Int)

    if nrow(df)!=1
        error("can fill in only a single dataframe row")
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
        error("can fill in only a single dataframe row")
    end
    dk = collect(keys(dict))
    for ik in dk
        dict[ik] = df[symbol(ik)][1]
    end

end

# dataframe to dict function
function df2dict(df::DataFrame)
  nm = names(df)
  snm = map(x->string(x),nm)
  out ={i => df[symbol(i)] for i in snm}
  return out
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

function fitMirror!(x::DataFrame,b::DataFrame)
    for i in 1:length(x)
        x[i] = fitMirror(convert(Float64,x[i][1]),convert(Float64,b[i,:lb][1]),convert(Float64,b[i,:ub][1]))
    end
end



# function checkbounds!(df::DataFrame,di::Dict)
#     if nrow(df) > 1
#         error("can only process a single row")
#     end
#     dfbounds = collectFields(di,1:length(di),true)
#     for c in names(df)
#         test1 = true 
#         test2 = true

#         while( test1 | test2)
#             # if above upper bound
#             if df[1,c] > dfbounds[2,c]
#                 df[1,c] = 2*dfbounds[2,c] - df[1,c]
#             else
#                 test1 = false
#             end
#             if df[1,c] < dfbounds[1,c]
#                 df[1,c] = 2*dfbounds[1,c] - df[1,c]
#             else
#                 test2 = false
#             end
#         end
#     end
#     return df
# end


            

