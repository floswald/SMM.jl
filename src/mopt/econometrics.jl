




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
# Compute Score of moments

Computes an approximation to the simulated score of moments at the optimal parameter value. This is the matrix of derivatives of moments w.r.t. parameters around the best parameter value.

"""
function score(MA::MAlgo)

	p = parameters_ID(MA.MChains)
	inf = infos(MA.MChains)
	best = findmax(inf[:evals])

	pnames = convert(Array{Symbol},sort(collect(keys(MA.m.params_to_sample))))

	p0= p[best[2],pnames]
	wts = mean((array(p[pnames]) .- array(p0)).^2,2)

	data = hcat(p,moments(MA.MChains))

	M  = Dict()
	lhs = sort(collect(keys(MA.m.moments)))

	for m in lhs
		fm = Formula(m, Expr(:call, :+, pnames...))
		s  = glm(fm, data, Normal(), IdentityLink(), wts = wts[:] )
		M[m] = coef(s)[2:end]
	end
	return M
end



"""
Computes standard errors of MCMC estimates
"""
function std(MA::MAlgo)

	Sd = score(MA)
	S = zeros(length(Sd),length(Sd[collect(keys(Sd))[1]]))
	row = 0
	for k in sort(collect(keys(Sd)))
		row +=1
		S[row,:] = Sd[k]
	end

	pnames = convert(Array{Symbol},sort(collect(keys(MA.m.params_to_sample))))

	# diagonal matrix of data sd's
	Sigma  = Diagonal( Float64[ MA.m.moments[ki][:weight] for ki in sort(collect(keys(Sd)))] )

	Omega = S * Sigma * S'

	# Covarance matrix of parameters from coldest chain
	chs = parameters_ID(MA.MChains)
	ch1 = chs[chs[:chain_id].==1,convert(Array{Symbol}, sort(collect(keys(MA.m.params_to_sample))))]
	J = cov(array(ch1))

	SE = J * Omega * J

	se = Dict()
	row=0
	for k in pnames
		row+=1
		se[k] = diag(SE)[row]
	end

	return se
end



