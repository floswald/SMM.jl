# We introduce two bojects to facilitate the interfacing
# with the objective function: Eval and Opts.

export Eval, start, finish, param, paramd, fill, dataMoment, dataMomentW, setMoment, setValue

type Eval

	value  ::Float64
	time   ::Float64
	params ::Dict
	moments::Dict
	dataMoments::Dict
	dataMomentsW::Dict
	status ::Int64
	options :: Dict

	function Eval(p::Dict,mom::DataFrame)
		this = new()
		this.value = -1
		this.time = time()
		this.status = -1
		this.dataMoments  = Dict{Symbol,Float64}()
		this.dataMomentsW = Dict{Symbol,Float64}()
		this.params       = Dict{Symbol,Float64}()
		this.moments      = Dict{Symbol,Float64}()

		for (i in 1:nrow(mom))
			kk = symbol(mom[i,:name])
			this.dataMoments[kk] = mom[i,:value]
			this.dataMomentsW[kk] = mom[i,:sd]
		end

		for (k in keys(p))
			kk = symbol(k)
			this.params[kk] = p[k]
		end

		return this
	end
end

function start(ev::Eval)
	ev.time = time()
end

function finish(ev::Eval)
	ev.time = ev.time - time()
end

param(ev::Eval,ll::Array{Symbol,1})    = [ ev.params[i] for i in ll]
param(ev::Eval,ll::Array{Any,1})       = [ ev.params[i] for i in ll]
param(ev::Eval)                        = param(ev,collect(keys(ev.params)))
param(ev::Eval,s::Symbol)              = param(ev,[s])
paramd(ev::Eval)                       = ev.params
dataMoment(ev::Eval,ll::Array{Symbol,1})  = [ ev.dataMoments[i] for i in ll]
dataMoment(ev::Eval)                   = dataMoment(ev,keys(ev.dataMoments))
dataMomentW(ev::Eval,ll::Array{Symbol,1}) = [ ev.dataMomentsW[i] for i in ll]
dataMomentW(ev::Eval)                  = dataMoment(ev,keys(ev.dataMomentsW))
dataMoment(ev::Eval,s::Symbol)         = dataMoment(ev,[s])
dataMomentW(ev::Eval,s::Symbol)        = dataMomentW(ev,[s])

# this allows to fill the values of a given structure
# with the values from ev
function fill(p::Any,ev::Eval)
	for k in keys(ev.params)
		setfield!(p,k,ev.params[k])
	end
end

function setValue(ev::Eval,value::Float64)
	ev.value = value
end

function setMoment(ev::Eval,k::Symbol,value::Float64)
	ev.moments[k] = value
end

function setMoment(ev::Eval,d::Dict)
	for k in keys(d)
		ev.moments[k] = d[k]
	end
end

# this assumes that colum :name has the names as strings
# and that column :value stores the value
function setMoment(ev::Eval,d::DataFrame)
	for i in 1:nrow(d)
		ev.moments[ symbol(d[i,:name]) ] = d[i,:value]
	end
end


