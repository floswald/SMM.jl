# MOpt

## Exported
---

#### Eval()
Eval Default Constructor: creates an empty `Eval` instance

**source:**
[MOpt/src/mopt/Eval.jl:29](file:///Users/florianoswald/.julia/v0.3/MOpt/src/mopt/Eval.jl)

---

#### Eval(mprob::MProb, p::Dict{K, V})
Eval Constructor for `MProb` and param dict:

* `MProb`: an `MProb` object
* `p`: Dict of parameters



**source:**
[MOpt/src/mopt/Eval.jl:82](file:///Users/florianoswald/.julia/v0.3/MOpt/src/mopt/Eval.jl)

---

#### Eval(p::Dict{K, V}, mom::DataFrame)
Eval Constructor for moment dataframe:

* `mom`: DataFrame with moments. Needs 3 columns named `name`, `value` and `weight`
* `p`: Dict of parameters



**source:**
[MOpt/src/mopt/Eval.jl:47](file:///Users/florianoswald/.julia/v0.3/MOpt/src/mopt/Eval.jl)

---

#### Eval(p::Dict{Symbol, Float64}, m::Dict{Symbol, (Float64, Float64)})
Eval Constructor for moment and param dicts:

* `m`: Dict of datamoments with weights
* `p`: Dict of parameters



**source:**
[MOpt/src/mopt/Eval.jl:130](file:///Users/florianoswald/.julia/v0.3/MOpt/src/mopt/Eval.jl)

---

#### Eval(p::Dict{Symbol, Float64}, m::Dict{Symbol, Float64})
Eval Constructor for moment and param dicts:

* `m`: Dict of datamoments without weights
* `p`: Dict of parameters



**source:**
[MOpt/src/mopt/Eval.jl:111](file:///Users/florianoswald/.julia/v0.3/MOpt/src/mopt/Eval.jl)

---

#### dataMoment(ev::Eval)
Obtain value(s) of selected moment(s) as an array

**source:**
[MOpt/src/mopt/Eval.jl:174](file:///Users/florianoswald/.julia/v0.3/MOpt/src/mopt/Eval.jl)

---

#### dataMoment(ev::Eval, ll::Array{Symbol, 1})
Obtain value(s) of selected moment(s) as an array

**source:**
[MOpt/src/mopt/Eval.jl:174](file:///Users/florianoswald/.julia/v0.3/MOpt/src/mopt/Eval.jl)

---

#### dataMoment(ev::Eval, s::Symbol)
Obtain value(s) of selected moment(s) as an array

**source:**
[MOpt/src/mopt/Eval.jl:174](file:///Users/florianoswald/.julia/v0.3/MOpt/src/mopt/Eval.jl)

---

#### dataMomentW(ev::Eval)
Obtain value(s) of selected moment weight(s) as an array

**source:**
[MOpt/src/mopt/Eval.jl:182](file:///Users/florianoswald/.julia/v0.3/MOpt/src/mopt/Eval.jl)

---

#### dataMomentW(ev::Eval, ll::Array{Symbol, 1})
Obtain value(s) of selected moment weight(s) as an array

**source:**
[MOpt/src/mopt/Eval.jl:182](file:///Users/florianoswald/.julia/v0.3/MOpt/src/mopt/Eval.jl)

---

#### dataMomentW(ev::Eval, s::Symbol)
Obtain value(s) of selected moment weight(s) as an array

**source:**
[MOpt/src/mopt/Eval.jl:182](file:///Users/florianoswald/.julia/v0.3/MOpt/src/mopt/Eval.jl)

---

#### fill(p, ev::Eval)
this allows to fill the values of a given structure with the values from Eval

**source:**
[MOpt/src/mopt/Eval.jl:190](file:///Users/florianoswald/.julia/v0.3/MOpt/src/mopt/Eval.jl)

---

#### param(ev::Eval)
Obtain value(s) of selected parameter(s) as an array

**source:**
[MOpt/src/mopt/Eval.jl:162](file:///Users/florianoswald/.julia/v0.3/MOpt/src/mopt/Eval.jl)

---

#### param(ev::Eval, ll::Array{Any, 1})
Obtain value(s) of selected parameter(s) as an array

**source:**
[MOpt/src/mopt/Eval.jl:162](file:///Users/florianoswald/.julia/v0.3/MOpt/src/mopt/Eval.jl)

---

#### param(ev::Eval, ll::Array{Symbol, 1})
Obtain value(s) of selected parameter(s) as an array

**source:**
[MOpt/src/mopt/Eval.jl:162](file:///Users/florianoswald/.julia/v0.3/MOpt/src/mopt/Eval.jl)

---

#### param(ev::Eval, s::Symbol)
Obtain value(s) of selected parameter(s) as an array

**source:**
[MOpt/src/mopt/Eval.jl:162](file:///Users/florianoswald/.julia/v0.3/MOpt/src/mopt/Eval.jl)

---

#### paramd(ev::Eval)
Obtain all paramter values as dict

**source:**
[MOpt/src/mopt/Eval.jl:166](file:///Users/florianoswald/.julia/v0.3/MOpt/src/mopt/Eval.jl)

---

#### Eval
`Eval` type for managing function evaluations 

## fields 

*`value`: current function value
*`time`: initial timing
*`params`: Dict of parameters
*`moments`: Dict of moments
*`dataMoments`: Dict of data moments
*`dataMomentsW`: Dict of weights for data moments
*`status`: Int error status
*`options`: Dict of options


**source:**
[MOpt/src/mopt/Eval.jl:16](file:///Users/florianoswald/.julia/v0.3/MOpt/src/mopt/Eval.jl)

---

#### MProb
This package implements several MCMC algorithms to optimize a non-differentiable objective function. The main application are **likelihood-free estimators**, which requires evaluating the objective at many regions. In general, this implements *Simulated Method of Moments*. The library is targeted to run MCMC on an SGE cluster, where each node is a chain.

**source:**
[MOpt/src/mopt/mprob.jl:11](file:///Users/florianoswald/.julia/v0.3/MOpt/src/mopt/mprob.jl)

## Internal
---

#### dataMomentWd(ev::Eval)
Obtain all moment weights as dict

**source:**
[MOpt/src/mopt/Eval.jl:186](file:///Users/florianoswald/.julia/v0.3/MOpt/src/mopt/Eval.jl)


