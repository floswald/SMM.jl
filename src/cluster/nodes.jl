
using MOpt

# data are generated from a bivariate normal
# with mu = [a,b] = [0,0]
# aim: 
# 1) sample [a',b'] from a space [-1,1] x [-1,1] and
# 2) find true [a,b] by computing distance(S([a',b']), S([a,b]))
#    and accepting/rejecting [a',b'] according to BGP
# 3) S([a,b]) returns a summary of features of the data

# initial value
p    = ["a" => 0.9 , "b" => -0.9]
# param bounds
pb   = [ "a" => [-1,1] , "b" => [-1,1] ]
# data moments
moms = MOpt.DataFrame(moment=["alpha","beta"],data_value=[0.0,0.0],data_sd=rand(2))

# define a minization problem
mprob = MProb(p,pb,MOpt.objfunc_norm2,moms)
