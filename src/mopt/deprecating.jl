
@deprecate(setMoment(ev::Eval,k::Symbol,value::Float64),setMoments!(ev::Eval,k::Symbol,value::Float64))
@deprecate(setMoment(ev::Eval,d::Dict),setMoments!(ev::Eval,d::Dict))
@deprecate(setMoment(ev::Eval,d::DataFrame),setMoments!(ev::Eval,d::DataFrame))
@deprecate(setValue(ev::Eval,x::Float64),setValue!(ev::Eval,x::Float64))