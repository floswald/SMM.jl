#'.. py:function:: Testobj
#'
#'   define a Test objective function
function Testobj(x::Dict,mom::DataFrame,whichmom::Array{String,1},vargs...)

	t0 = time()

    if length(vargs) > 0
        if get(vargs[1],"printlevel",0) > 0
            Base.info("in Test objective function")
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
	ret = Dict("value" => v, "params" => deepcopy(x), "time" => t0, "status" => status, "moments" => mout)
	return ret

end


export Testobj2,Testobj3,objfunc_norm

# dummy objective function
# this does not have a well defined minimum
# so will not work for estimation
function Testobj2(ev::Eval)
    start(ev)
    # info("in Test objective function Testobj2")

    val = 0.0
    for (k,v) in dataMomentd(ev)
        setMoment(ev,k,v+2.2)
        val += (2.2)^2
    end

	ev.status = 1
    setValue(ev,val)
    return ev
end


function Testobj_fails(ev::Eval)
	
	# this function returns with an exception
	open("/no/data/here")

end

#'.. py:function:: objfunc_norm
#'
#'   define a Test objective function
function objfunc_norm(ev::Eval)
    
	start(ev)
	# info("in Test objective function objfunc_norm")

	# extract parameters    
	mu  = convert(Array{Float64,1},param(ev)) # returns entire parameter vector 
	# use paramd(ev) to get as a dict.

	# compute simulated moments
	ns = 5000
	sigma           = convert(Matrix,Diagonal([1.0,1.0]))
	randMultiNormal = MOpt.MvNormal(mu,sigma) 
	simM            = mean(rand(randMultiNormal,ns),2)
	simMoments = Dict(:mu1 => simM[1], :mu2 => simM[2])


	# get data mometns
	# same thing here, and use dataMomentsWeights for sd
	# second argument can be optional
	# get objective value: (data[i] - model[i]) / weight[i]
	v = Dict{Symbol,Float64}()
	for (k,mom) in dataMomentd(ev)
		if haskey(dataMomentWd(ev),k)
			v[k] = ((simMoments[k] .- mom) ./ dataMomentW(ev,k)) .^2
		else
			v[k] = ((simMoments[k] .- mom) ) .^2
		end
	end
	setValue(ev, mean(collect(values(v))))
	# value = data - model
	# setValue(ev, mean((simMoments - trueMoments).^2) )

	# also return the moments
	setMoment(ev, simMoments)
	# mdf = DataFrame(name=["m1","m2"],value=simMoments[:])
	# setMoment(ev, mdf)

	ev.status = 1

	# finish and return
	finish(ev)

    return ev
end



#'.. py:function:: banana
#'
#'   define a Test objective function
function banana(ev::Eval)

    start(ev)
	p = paramd(ev)
    model = 100 .* (p["b"] - p["a"].^2 ).^2 .+ (1.-p["a"])^2
    data  = 0.0

    setValue(ev,model)
    for (k,v) in dataMomentd(ev)
        setMoment(ev,k,v+2.2)
    end
    finish(ev)
    return ev

end