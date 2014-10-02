module TestBGPChain

using FactCheck, DataFrames, MOpt


# TESTING Chains
# ==============
p    = [:a => 3.1 , :b => 4.9]
moms = DataFrame(name=[:alpha,:beta,:gamma],value=[0.8,0.7,0.5],sd=rand(3))

p2    = ["a" => 3.1 , "b" => 4.9]
moms2 = DataFrame(name=[:alpha,:beta,:gamma],value=[0.8,0.7,0.5],sd=[0.1,0.2,0.3])

type MyP
	a :: Float64 
	b :: Float64 

	function MyP()
		return  new(0.0,0.0)
	end
end


facts("Testing Eval object") do

	ev= Eval(p,moms)
	ev2= Eval(p2,moms2)	


	context("testing param") do
		@fact param(ev,:a) => [3.1]
		@fact param(ev,[:a,:b]) => [3.1,4.9]

		@fact param(ev2,:a) => [3.1]
		@fact param(ev2,[:a,:b]) => [3.1,4.9]
		@fact paramd(ev) => [:a => 3.1 , :b => 4.9]

		myp = MyP()
		fill(myp,ev)
		@fact myp.a => 3.1
	end

	context("testing moments") do
		@fact dataMoment(ev2,:alpha) => [0.8]
		@fact dataMomentW(ev2,:alpha) => [0.1]	

		setMoment(ev,:alpha,0.78)	
		setMoment(ev2,DataFrame(name=["alpha"] , value = [0.78]))	

		@fact ev.moments[:alpha] => 0.78
		@fact ev2.moments[:alpha] => 0.78

		setValue(ev,4.2)
		@fact ev.value => 4.2
	end

end



end # module 