



# minimize the rosenbrock function
# y = 100 * (b - a * a)^2 + (1 - a)^2
function banana(x::Dict,mom::Dict,whichmom::Array{ASCIIString,1})
	v = 100 * (x["b"] - x["a"]^2 ) ^2 + (1-x["a"])^2
	# we just want to find the lowest value - no moments involved.
	ret = ["value" => v, "params" => x, "time" => 1.0, "status" => 1, "moments" => mom]
	return ret
end


p    = ["a" => 3.1 , "b" => 4.9]
pb   = [ "a" => [-100,100] , "b" => [-100,100] ]
moms = [
	"alpha" => [ 0.8 , 0.02 ],
	"beta"  => [ 0.8 , 0.02 ],
	"gamma" => [ 0.8 , 0.02 ]
]



mprob = Mopt.MProb(p,pb,banana,moms)
opts =["N"=>20,"shock_var"=>1.0,"mode"=>"serial","maxiter"=>100,"path"=>".","maxtemp"=>100,"min_shock_sd"=>1.0,"max_shock_sd"=>15.0,"past_iterations"=>30,"min_jumptol"=>0.1,"max_jumptol"=>10.0] 
MA = Mopt.MAlgoBGP(mprob,opts)

Mopt.runMopt(MA)


x=Mopt.alls(MA.MChains)
x[x[:chain_id] .== 1,:]
println(Mopt.alls(MA.MChains))




