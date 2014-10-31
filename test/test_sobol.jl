module TestSobol

	using FactCheck,DataFrames,Lazy,MOpt

	# initial value
	pb    = ["m1" => [0.0,-2,2] , "m2" => [0.0,-2,2] ] 
	moms = DataFrame(name=["m1","m2"],value=[0.0,0.0],weight=rand(2))
	mprob = @> MOpt.MProb() MOpt.addSampledParam!(pb) MOpt.addMoment(moms) MOpt.addEvalFunc(MOpt.objfunc_norm);

	# compute the slices
	sl = MOpt.sobolsearch(mprob,500);

	# find the min value

	best=sl[1]
	for e in sl
		if e.value < best.value
			best=e
		end
	end

	print(best.params)


	sl = MOpt.sobolWeightedSearch(mprob,500);

	using PyPlot

	plot(sl[:X][:,1],sl[:X][:,2],"o")


	# subplot(231)
	# r = MOpt.get(sl, :m1 , :m1); PyPlot.plot(r[:x],r[:y],".")
	# subplot(232)
	# r = MOpt.get(sl, :m1 , :m2); PyPlot.plot(r[:x],r[:y],".")
	# subplot(233)
	# r = MOpt.get(sl, :m1 , :value); PyPlot.plot(r[:x],r[:y],".")
	# subplot(234)
	# r = MOpt.get(sl, :m2 , :m1); PyPlot.plot(r[:x],r[:y],".")
	# subplot(235)
	# r = MOpt.get(sl, :m2 , :m2); PyPlot.plot(r[:x],r[:y],".")
	# subplot(236)
	# r = MOpt.get(sl, :m2 , :value); PyPlot.plot(r[:x],r[:y],".")


	# generate a sobol sequence
	using Sobol
	sq = SobolSeq(2)
	S = [ next(sq) for k in 1:10000]
	S1 = [ x[1] for x in S]
	S2 = [ x[2] for x in S]

	# define a simple objective function
	mu = [1.0,1.0]
	C  = [0.5 0.3; 0.1 0.8]
	R = MvNormal(mu, C)
	function f(x)
	 x1 = quantile(Normal(),x[1])
	 x2 = quantile(Normal(),x[2]) 
	 return pdf(R,[x1,x2])
	end

	V = map(f,S)

	PyPlot.scatter(S1,S2,c=V)

	# how to decide to add new points?

	# compute 2 things at each point: 
	# - a notion of density
	# - a notion of level

	# evaluate the first 100 points
	E = V*0
	E[1:100] = 1

	# given the evaluated set, compute the local level and densities
	N = sum(E)^(1/dim)
    h = 1/N # distance between points

    # for each point, we compute a local average of the value
    K  = length(S)
    S1e = S1[E .== 1] 
    S2e = S2[E .== 1] 
    Se  = S[E .== 1] 
    Ve  = V[E .== 1] 
    Sa = [S1 S2]


    # take 100 evaluated points and compute V and D 
    # for them and average.

    # add a point to Sobol list, check if it should be evaluated
    # randomly pick an unevaluated point from the past check if it should be evaluated

    # pick a new evaluation point, based on previously computed.
    sq   = SobolSeq(2)
    P0 = [ [ :p => next(sq) , :value => 0] for k in 1:10]

    # evaluate 10 points:
    P1 = map(Pt0) do p
    	p[:value] = f(p[:p])
    	p
    end

    # get an approx distance
	N = length(P1)^(1/dim)
    h = 1/N # distance between points


    # 1 find closest evalued points
    function findClose(s,P)
    	best_dist = Inf
    	sbest = P[1]
	    for s2 in P
	    	d = sum( (s.-s2[:p]).^2)
	    	if (  d > 0 ) & ( d < best_dist )
	    		sbest = s2
	    		best_dist = d
	    	end
	    end
		return (sbest,best_dist)
	end




    d = [ exp(  - sum( ( S[i] .- ss ).^2 )/h^2 ) for ss in Se ]




    for i in 1:K
    	w = [ exp(  - sum( ( S[i] .- ss ).^2 )/h^2 ) for ss in Se ]
    	V_hat = sum ( w .* Ve) ./ sum( w )
    	w = [ exp(  - sum( ( S[i] .- ss ).^2 )/h^2 ) for ss in S ]
    	D_hat = sum( w .* E) ./ sum (w)
    end


    # compute the density
 









end # module 






