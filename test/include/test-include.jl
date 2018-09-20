
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
    bundle = 2
    chain = MomentOpt.BGPChain(id,n,m=mprob,sig=sig,upd=upd,upd_by=upd_by,smpl_iters=ite,upd_bundle=bundle)
    (chain, id, n, mprob, sig, sig2, upd, upd_by, ite,bundle)
end

function test_chain2()
    pb   = Dict( "a" => [0.3; -1;1] , "b" => [-0.9;-2;2] , "c" => [-0.9;-2;2], "d" => [-0.9;-2;2], "e" => [-0.9;-2;2])
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
    bundle = 2
    chain = MomentOpt.BGPChain(id,n,m=mprob,sig=sig,upd=upd,upd_by=upd_by,smpl_iters=ite,upd_bundle=bundle)
    (chain, id, n, mprob, sig, sig2, upd, upd_by, ite,bundle)
end

type MyP
    a :: Float64 
    b :: Float64 

    function MyP()
        return  new(0.0,0.0)
    end
end