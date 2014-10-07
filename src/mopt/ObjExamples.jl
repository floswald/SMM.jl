#'.. py:function:: Testobj
#'
#'   define a Test objective function
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


export Testobj2,Testobj3,objfunc_norm

function Testobj2(ev::Eval)

    start(ev)
    info("in Test objective function")

    val = 0.0
    for (k,v) in dataMoment(ev)
        setMoment(ev,k,v+2.2)
        val += (2.2)^2
    end

    setValue(ev,val)
    return ev
end

#'.. py:function:: objfunc_norm2
#'
#'   define a Test objective function
function objfunc_norm(ev::Eval)
    
	start(ev)
	info("in Test objective function")

	# extract parameters    
	mu  = convert(Array{Float64,1},param(ev,[:m1,:m2])) # returns values as a vector, use param to get a Dict

	# compute simulated moments
	ns = 5000
	sigma           = convert(Matrix,Diagonal([1.0,1.0]))
	randMultiNormal = MOpt.MvNormal(mu,sigma) 
	simMoments      = mean(rand(randMultiNormal,ns),2)

	# get data mometns
	trueMoments = dataMoment(ev,[:m1,:m2]) 
	# same thing here, and use dataMomentsWeights for sd
	# second argument can be optional

	# value = data - model
	setValue(ev, mean((simMoments - trueMoments).^2) )

	# also return the moments
	setMoment(ev, {:m1 => simMoments[1], :m2 => simMoments[2]})
	# we would also have a setter that takes a DataFrame

	ev.status = 1

	# finish and return
	finish(ev)

    return ev
end

#'.. py:function:: banana
#'
#'   define a Test objective function
function banana(x::Dict,mom::DataFrame,whichmom::Array{ASCIIString,1})

    model = 100 .* (x["b"] - x["a"].^2 ).^2 .+ (1.-x["a"])^2
    data  = 0.0

    value = mean((data .- model).^2)

    momout = DataFrame(alpha = model)



    # we just want to find the lowest value - no moments involved.
    ret = ["value" => value, "params" => x, "time" => 1.0, "status" => 1, "moments" => momout]
    return ret
end