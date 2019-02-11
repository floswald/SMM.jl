




"""
load algo from file and compute simple standard errors as in Lise, Meghir, Robin (2016)
Because of lack of numerical precision in numerical approx of gradient, function [`score`](@ref MOpt.score), this reports simply the standard deviation of the coldest chain.
"""
function get_simple_std(f::String)
	x = readAlgoBGP(f)
	nms = setdiff(names(x["params"]),[:chain_id,:iter])
	m = colwise(mean,@where(x["params"],:chain_id.==1)[nms])
	s = colwise(std,@where(x["params"],:chain_id.==1)[nms])
	m = Float64[m[im][1] for im in 1:length(m)]
	s = Float64[s[im][1] for im in 1:length(s)]
	out = DataFrame(param=string.(nms),estimate=m,sd=s)
end


"""
	FD_gradient(m::MProb,p::Dict;step_perc=0.005)

Get the gradient of the moment function wrt to some parameter vector via finite difference approximation. 
The output is a (k,n) matrix, where ``k`` is the number of `m.params_to_sample` and where ``m`` is the number of moments.

The default step size is 1% of the parameter range.
"""
function FD_gradient(m::MProb,p::Union{Dict,OrderedDict};step_perc=0.01)

	# get g(p)
    ev = evaluateObjective(m,p)
    mnames = collect(keys(m.moments))
    smm = filter((x,y)->in(x,mnames),ev.simMoments)
	gp = collect(values(smm))
	D = zeros(length(p),length(gp))

	# optimal step size depends on range of param bounds
	rs = range_length(m)

	# compute each partial derivative
	rows = pmap( [(k,v) for (k,v) in p ] ) do ip 
		k = ip[1]
		v = ip[2]
		h = rs[k] * step_perc
		pp = deepcopy(p)
		pp[k] = v + h 
		# println("changing $k from $v to $(pp[k]) by step $h")
		xx = evaluateObjective(m,pp)
		smm = collect(values(filter((x,y)->in(x,mnames),xx.simMoments)))
		Dict(:p => k, :smm => (smm .- gp) / h)
	end
	d = Dict()
	for e in rows
       d[e[:p]] = e[:smm]
    end
	row = 0
	for (k,v) in d
		row += 1
		D[row,:] = v
	end

	return D

end

function get_stdErrors(m::MProb,p::Union{Dict,OrderedDict};reps=300)

	# compute "Data" var-cov matrix Σ by generating H samples of simulated data using p 
	Σ = getSigma(m,p,reps)

	# compute score of moment function
	J = FD_gradient(m,p)
	println("size of gradient J = $(size(J))")

	# put all together to get standard errors
	# first get weighting matrix 
	W = Diagonal([v[:weight] for (k,v) in m.moments])
	println("size of Weight W = $(size(W))")
	SE = pinv(J *  W * J') * (J * W * Σ * W * J') * pinv(J * W * J') 
	return Dict(zip(collect(keys(p)),sqrt.(diag(SE))))

end



