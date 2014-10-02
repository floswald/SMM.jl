# We introduce two bojects to facilitate the interfacing
# with the objective function: Eval and Opts.

type Eval

	value  ::Float64
	time   ::Float64
	params ::Dict
	moments::Dict
	dataMoments::Dict
	status ::Int64

	function Eval(p::Dict,mom::DataFrame)
		this = new()
		this.value = -1
		this.time = time()
		this.params = deepcopy(p)
		this.status = -1
		this.moments = mom
	end
end

