
# define a Test objective function
function Testobj(x::Dict,mom::DataFrame,whichmom::Array{ASCIIString,1},vargs...)

	t0 = time()

    if length(vargs) > 0
        if get(vargs[1],"printlevel",0) > 0
            info("in Test objective function")
        end
    end

    mm = deepcopy(mom)
    nm0 = names(mm)
    DataFrames.insert_single_column!(mm,zeros(nrow(mm)),ncol(mm)+1)
    names!(mm,[nm0,:model_value])

    for ir in eachrow(mm)
        ir[:model_value] = ir[:data_value] + 2.2
    end

	# output all moments
    mout = transpose(mm[[:moment,:model_value]],1)

	# subset mm to required moments to compute function value
	mm = mm[findin(mm[:moment],whichmom),:]

	# compute distance
	v = sum((mm[:data_value] - mm[:model_value]).^2)

	# status
	status = 1

	# time out
	t0 = time() - t0

	# return a dict
	ret = ["value" => v, "params" => deepcopy(x), "time" => t0, "status" => status, "moments" => mout]
	return ret

end

 # ret = ["value" => 1.1, "params" => ["a"=>1.1,"b"=>12.1], "time" => 0, "status" => 1, "moments" => ["alpha"=>1.1,"beta"=>12.1,"gamma"=>12.1] ]


function objfunc_norm2(p::Dict,mom::DataFrame,whichmom::Array{ASCIIString,1},vargs...)
    
    # println("I am worker number: $(myid())")

    t0 = time()
    if length(vargs) > 0
        if get(vargs[1],"printlevel",0) > 0
            info("in Test objective function")
        end
    end
    
    sigma = convert(Matrix,Diagonal([1.0,1.0]))
    ns = 5000

    # compute simulated moments
    mu = [p["a"],p["b"]]
    MVN = MOpt.MvNormal(mu,sigma) 

    # get data mometns
    muD = array(transpose(mom[:,[1,2]],1))'

    # simulate model moments 
    moments = mean(rand(MVN,ns),2)
        
    # value = data - model
    value = mean((muD - moments).^2)

    momout = DataFrame(alpha = moments[1],beta = moments[2])

    t0 = time() - t0

        ret = ["value"=>value, "params" =>p, "time" =>t0 , "status" => 1, "moments" => momout]
    return ret
end



function banana(x::Dict,mom::DataFrame,whichmom::Array{ASCIIString,1})

    model = 100 .* (x["b"] - x["a"].^2 ).^2 .+ (1.-x["a"])^2
    data  = 0.0

    value = mean((data .- model).^2)

    momout = DataFrame(alpha = model)



    # we just want to find the lowest value - no moments involved.
    ret = ["value" => value, "params" => x, "time" => 1.0, "status" => 1, "moments" => momout]
    return ret
end

# transpose a 2-column dataframe to a one-row dataframe, so that 
# col x become new column names
function transpose(x::DataFrame,newNames::Int)
    if ncol(x) != 2
        throw(ArgumentError("x must have 2 columns"))
    end
    if !in(newNames,[1,2])
        throw(ArgumentError("newNames indexes col with new names: either 1 or 2"))
    end
    newrows = setdiff([1,2],newNames)[1]
    # make new column names out of col index newNames
    z = DataFrame(Float64,1,nrow(x))
    names!(z,Symbol[y for y in x[:,newNames]])
    for i in 1:nrow(x)
        z[1,i] = x[i,newrows]
    end
    return z
end



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
    dk = names(df)
    for ik in dk
        dict[string(ik)] = df[ik][1]
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


            

