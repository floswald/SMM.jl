# smart grid search

module AGrids

type Vertex

	x      ::Array{Float64,1}
	value  ::Float64
	closest_d :: Float64
	closest_v :: Float64

	function Vertex(x)
		this   = new()
		this.x = x
		this.value = Inf
		this.closest_d = Inf
		this.closest_v = Inf
		return this
	end
end

type AGrid
	points     :: Array{Vertex,1}

	function AGrid()
		this = new()
		this.points = Vertex[] 
		return this
	end
end

function findClose(ag::AGrid,x)
	best_dist = Inf
	sbest = P[1]
	i = 0
    for s2 in ag.points
    	i = i+1
    	d = sum( (s.- s2.x).^2)
    	if (  d > 0 ) & ( d < best_dist )
    		sbest = s2
    		best_dist = d
    	end
    end
	return (i,sbest,best_dist)
end

# adds a new evaluation to the grid
# and updates the closest neighbor for each guy
function add!(ag::AGrid,x,v)
	closest_d = Inf
	closest_v = Inf

	for s2 in ag.points
	  d = sum( ( x.- s2.x ).^2 )
	  
	  # update closest for the old
	  if d < s2.closest_d
	  	s2.closest_d = d
	  	s2.closest_v = v
	  end

	  # update clostest for current
	  if d < closest_d
	  	closest_d = d
	  	closest_v = s2.value
	  end

	end

	vv = Vertex(x)
	vv.value = v
	vv.closest_d = closest_d
	vv.closest_v = closest_v
	push!(ag.points,vv)
end

function energy(v::Vertex)
	(p.closest_d).^2 - p.closest_v
end


# compute the median on a combination 
# of discepency and level
function cutoff(ag::AGrid)
	return median([ energy(p) for p in ag.points])
end

# decides if value should be evaluated
function shouldEval(ag::AGrid, x::Array{Float64,1})
	(d,v) = 




end


end