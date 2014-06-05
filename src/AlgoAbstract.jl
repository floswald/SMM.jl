

type ChainAlgoAbstract

end

function computeNewGuess( algo::ChainAlgoAbstract  )

end


type ChainAlgoRandom <: ChainAlgoAbstract
  i :: int

  function ChainAlgoRandom ()
    return new(0)
  end
end

function computeNewGuess( algo::ChainRandom  )
  algo.i = algo.i +1
end


