

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

function BGP_approx_posterior(theta::Float64)
    eps = 0.025
    N = Normal()
    e1 = eps - theta
    e2 = -eps - theta
    0.45 * ( cdf(N,e1) + cdf(N,10*e1) - cdf(N,e2) - cdf(N,10*e2)) + 0.1 * (cdf(N,e1 + 5) - cdf(N,e2 + 5))
end

function objfunc_BGP(ev::Eval)
    start(ev)
    p = paramd(ev)  # param vector as dict
    y = BGP_approx_posterior(p[:theta])
    # no moments to set...
    setValue(ev,-y)  # want to minimize
    ev.status=1
    finish(ev)
    return(ev)
end


"""
    objfunc_norm(ev::Eval)

Test objective function. This is a bivariate normal distribution that simumlates data from the parameters in `ev`. From this simulated data, sample means are computed, which should be close to the empirical moments on `ev`. The aim is to minimize this function.
"""
function objfunc_norm(ev::Eval)
    
	start(ev)
	# info("in Test objective function objfunc_norm")

	# extract parameters    
    # mu  = convert(Array{Float64,1},param(ev)) # returns entire parameter vector 
	mu  = collect(values(ev.params))
	# use paramd(ev) to get as a dict.

	# compute simulated moments
	ns = 10000
	sigma           = [1.0;1.0]
	randMultiNormal = MomentOpt.MvNormal(mu,MomentOpt.PDiagMat(sigma)) 
	simM            = mean(rand(randMultiNormal,ns),dims = 2)
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



function banana(ev::Eval)

    start(ev)
	p = paramd(ev)
    model = 100 .* (p["b"] - p["a"].^2 ).^2 .+ (1 .- p["a"])^2
    data  = 0.0

    setValue(ev,model)
    for (k,v) in dataMomentd(ev)
        setMoment(ev,k,v+2.2)
    end
    finish(ev)
    return ev

end