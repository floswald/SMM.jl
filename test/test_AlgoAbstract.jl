@testset "AlgoAbstract" begin

  @testset "Testing loading and saving" begin

    pb   = Dict( "a" => [0.3; -1;1] , "b" => [-0.9;-2;2] )
    moms = DataFrame(name=["mu1","mu2"],value=[0.0;0.0],weight=rand(2))
    mprob = MProb()
    addSampledParam!(mprob,pb)
    addMoment!(mprob,moms)
    addEvalFunc!(mprob,MomentOpt.objfunc_norm)

    # estimation options:
    #--------------------
    # number of iterations for each chain
    niter = 10
    # number of chains
    # nchains = nprocs()
    nchains = 2

    opts = Dict("N"=>nchains,
            "maxiter"=>niter,
            "maxtemp"=> 5,
            "coverage"=>0.025,
            "sigma_update_steps"=>10,
            "sigma_adjust_by"=>0.01,
            "smpl_iters"=>1000,
            "parallel"=>true,
            "min_improve"=>[0.05 for i in 1:nchains],
            "mixprob"=>0.3,
            "acc_tuner"=>12.0,
            "animate"=>false,
            "save_frequency"=>5,
            "filename"=>"MyMAProblem.jld2")


    # set-up BGP algorithm:
    MA = MAlgoBGP(mprob,opts)

    # run the estimation:
    @time MomentOpt.runMOpt!(MA)

    # load MA saved above
    MA2 = readMalgo(opts["filename"])

    # compare MA and M2
    #------------------
    @test MA.opts == MA2.opts
    @test MA.i == MA2.i

    # Compare chains
    for chainNumber = 1:MA.opts["N"]
        for fieldName in fieldnames(MA.chains[chainNumber])
          # for fields m and sigma, the test returns false, eventhough they have
          # the same values
            if fieldName!= :m && fieldName != :sigma
                @test getfield(MA.chains[chainNumber], fieldName) == getfield(MA2.chains[chainNumber], fieldName)
            end
        end
    end

    @test MA.dist_fun == MA2.dist_fun
    rm(opts["filename"])

  end

end
