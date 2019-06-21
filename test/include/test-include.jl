
function test_chain()
    pb   = Dict( "a" => [0.3; -1;1] , "b" => [-0.9;-2;2] )
    moms = DataFrame(name=["mu1","mu2"],value=[0.0;0.0],weight=rand(2))
    mprob = MProb() 
    addSampledParam!(mprob,pb) 
    addMoment!(mprob,moms) 
    addEvalFunc!(mprob,MomentOpt.objfunc_norm)
    id = 180
    n = 23
    sig = rand()
    sig2 = rand()
    upd = 5
    upd_by = rand()
    ite = 1000
    b_size = 2
    chain = MomentOpt.BGPChain(id,n,m=mprob,sig=sig,upd=upd,upd_by=upd_by,smpl_iters=ite,batch_size=b_size)
    (chain, id, n, mprob, sig, sig2, upd, upd_by, ite,b_size)
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
    sig = rand()
    sig2 = rand()
    upd = 5
    upd_by = rand()
    ite = 1000
    b_size = 2
    chain = MomentOpt.BGPChain(id,n,m=mprob,sig=sig,upd=upd,upd_by=upd_by,smpl_iters=ite,batch_size=b_size)
    (chain, id, n, mprob, sig, sig2, upd, upd_by, ite,b_size)
end

function test_chain3()
    pb = OrderedDict()
    pb["p1"] = [0.2,0,1]
    pb["p2"] = [0.002,0,0.02]
    pb["p3"] = [-0.3,-2,2]
    pb["p4"] = [-0.4,-2,2]
    pb["p5"] = [0.008,0,0.02]
    pb["p6"] = [0.4,-2,2]
    pb["p7"] = [0.2,0,1]
    pb["p8"] = [0.002,0,0.02]
    pb["p9"] = [-0.3,-2,2]
    pb["p10"] = [-0.4,-2,2]
    pb["p11"] = [0.008,0,0.02]
    pb["p12"] = [0.4,-2,2]
    pb["p13"] = [0.2,0,1]
    pb["p14"] = [0.002,0,0.02]
    pb["p15"] = [-0.3,-2,2]
    pb["p16"] = [-0.4,-2,2]
    pb["p17"] = [0.008,0,0.02]
    pb["p18"] = [0.4,-2,2]
    moms = DataFrame(name=["mu1",
                           "mu2",
                           "mu3",
                           "mu4",
                           "mu5",
                           "mu6",
                           "mu7",
                           "mu8",
                           "mu9",
                           "mu10",
                           "mu11",
                           "mu12",
                           "mu13",
                           "mu14",
                           "mu15",
                           "mu16",
                           "mu17",
                           "mu18"],
        value=[-1.0,0.0,0.5,-0.5,0.008,-0.7,-1.0,0.0,0.5,-0.5,0.008,-0.7,-1.0,0.0,0.5,-0.5,0.008,-0.7],weight=[-1.0,0.0,0.5,-0.5,0.008,-0.7,-1.0,0.0,0.5,-0.5,0.008,-0.7,-1.0,0.0,0.5,-0.5,0.008,-0.7])
    mprob = MProb() 
    addSampledParam!(mprob,pb) 
    addMoment!(mprob,moms) 
    addEvalFunc!(mprob,MomentOpt.objfunc_norm)
    id = 180
    n = 100
    sig = rand()
    sig2 = rand()
    upd = 5
    upd_by = rand()
    ite = 1000
    b_size = 1
    chain = MomentOpt.BGPChain(id,n,m=mprob,sig=sig,upd=upd,upd_by=upd_by,smpl_iters=ite,batch_size=b_size)
    (chain, id, n, mprob, sig, sig2, upd, upd_by, ite,b_size)
end

mutable struct MyP
    a :: Float64 
    b :: Float64 

    function MyP()
        return  new(0.0,0.0)
    end
end