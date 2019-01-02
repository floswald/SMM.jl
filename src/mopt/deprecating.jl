function setMoment(ev::Eval,k::Symbol,value::Float64)
	ev.simMoments[k] = value
end

function setMoment(ev::Eval,d::Dict)
	for k in keys(d)
		ev.simMoments[k] = d[k]
	end
end

# this assumes that colum :name has the names as strings
# and that column :value stores the value
function setMoment(ev::Eval,d::DataFrame)
	for i in 1:nrow(d)
		ev.simMoments[ Symbol(d[i,:name]) ] = d[i,:value]
	end
end

@deprecate(setMoment(ev::Eval,k::Symbol,value::Float64),setMoments!(ev::Eval,k::Symbol,value::Float64))
@deprecate(setMoment(ev::Eval,d::Dict),setMoments!(ev::Eval,d::Dict))
@deprecate(setMoment(ev::Eval,d::DataFrame),setMoments!(ev::Eval,d::DataFrame))
@deprecate(setValue(ev::Eval,x::Float64),setValue!(ev::Eval,x::Float64))