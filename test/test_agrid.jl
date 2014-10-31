include("src/mopt/agrid.jl")

using Distributions
using Sobol

sq   = SobolSeq(2)
P0 = [ [ :p => next(sq) , :value => 0] for k in 1:10]

# define a simple objective function
mu = [1.0,1.0]
C  = [0.5 0.3; 0.1 0.8]
R = MvNormal(mu, C)
function f(x)
	x1 = quantile(Normal(),x[1])
	x2 = quantile(Normal(),x[2]) 
	return pdf(R,[x1,x2])
end


ag = AGrids.AGrid()

x = next(sq)
v = f(x)
AGrids.add!(ag,x,v)