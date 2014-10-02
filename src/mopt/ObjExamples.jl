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


#'.. py:function:: objfunc_norm2
#'
#'   define a Test objective function
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


#'.. py:function:: objfunc_norm2
#'
#'   define a Test objective function
function objfunc_norm3(ev::Eval)
    
	start(ev)
	info("in Test objective function")

	# extract parameters    
	mu  = param(ev,[:a,:b]) # returns values as a vector, use param to get a Dict

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
	setMoment(ev, {:m1 => moments[1], :m2 => moments[2]})
	# we would also have a setter that takes a DataFrame

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