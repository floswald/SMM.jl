
function test_chain()
    pb   = Dict( "a" => [0.3; -1;1] , "b" => [-0.9;-2;2] )
    moms = DataFrame(name=["mu1","mu2"],value=[0.0;0.0],weight=rand(2))
    mprob = MProb() 
    addSampledParam!(mprob,pb) 
    addMoment!(mprob,moms) 
    addEvalFunc!(mprob,MomentOpt.objfunc_norm)
    id = 180
    n = 23
    sig = rand(length(pb))
    sig2 = rand(length(pb))
    upd = 5
    upd_by = rand()
    ite = 1000
    chain = MomentOpt.BGPChain(id,n,mprob,sig,upd,upd_by,ite)
    (chain, id, n, mprob, sig, sig2, upd, upd_by, ite)
end


mutable struct MyP
    a :: Float64 
    b :: Float64 

    function MyP()
        return  new(0.0,0.0)
    end
end