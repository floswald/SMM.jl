

Writing an objective function
-----------------------------

We show here how to write your own objective function. The signature is:

    obj(x::Dict,mom::DataFrame,whichmom::Array{ASCIIString,1},vargs...)

where the first agrument is a dictionary of values for the parameters. The second is a dataFrame 
that contains information about the different moments to match. This will include the value of the moments
together with their standard error.


The objective function should return another Dictionary that includes:

 - `value`: the value of the objective function
 - `params`: the list of parameters with their value, this should just be a deep copy of x
 - `time`: the amount of time spent in the objective function
 - `moments`: a DataFrame with the evaluated moments

here is a full example:

.. code-block:: julia

  function Testobj(x::Dict,mom::DataFrame,whichmom::Array{ASCIIString,1},vargs...)

    t0 = time()

    if length(vargs) > 0
        if get(vargs[1],"printlevel",0) > 0
            info("in Test objective function")
        end
    end

    mm = deepcopy(mom)
    nm0 = names(mm)
    DataFrames.insert_single_column!(mm,zeros(nrow(mm)),ncol(mm)+1)
    names!(mm,[nm0,:model_value])

    for ir in eachrow(mm)
        ir[:model_value] = ir[:data_value] + 2.2
    end

    # output all moments
    mout = transpose(mm[[:moment,:model_value]],1)

    # subset mm to required moments to compute function value
    mm = mm[findin(mm[:moment],whichmom),:]

    # compute distance
    v = sum((mm[:data_value] - mm[:model_value]).^2)

    # status
    status = 1

    # time out
    t0 = time() - t0

    # return a dict
    ret = ["value" => v, "params" => deepcopy(x), "time" => t0, "status" => status, "moments" => mout]
    return ret

  end




