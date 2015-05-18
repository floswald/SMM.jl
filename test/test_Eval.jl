module TestBGPChain

using FactCheck, DataFrames, MOpt, HDF5



# TESTING Eval 
# ============
p    = [:a => 3.1 , :b => 4.9]
moms = DataFrame(name=[:alpha,:beta,:gamma],value=[0.8,0.7,0.5],weight=rand(3))

p2    = ["a" => 3.1 , "b" => 4.9]
moms2 = DataFrame(name=[:alpha,:beta,:gamma],value=[0.8,0.7,0.5],weight=[0.1,0.2,0.3])

type MyP
	a :: Float64 
	b :: Float64 

	function MyP()
		return  new(0.0,0.0)
	end
end


facts("Testing Eval object") do

	ev = Eval(p,moms)
	setMoment(ev,moms)
	ev.status=1
	ev2= Eval(p2,moms2)	
	setMoment(ev2,moms2)
	ev2.status=1

	context("testing param") do

		@fact param(ev,:a) => 3.1
		@fact param(ev,[:a,:b]) => [3.1,4.9]

		@fact param(ev2,:a) => 3.1
		@fact param(ev2,[:a,:b]) => [3.1,4.9]
		@fact paramd(ev) => [:a => 3.1 , :b => 4.9]

		myp = MyP()
		fill(myp,ev)
		@fact myp.a => 3.1
	end

	context("testing moments") do
		
		@fact dataMoment(ev2,:alpha) => 0.8
		@fact dataMomentW(ev2,:alpha) => 0.1	

		setMoment(ev,:alpha,0.78)	
		setMoment(ev2,DataFrame(name="alpha" , value = 0.78))	

		@fact ev.simMoments[:alpha] => 0.78
		@fact ev2.simMoments[:alpha] => 0.78

		setMoment(ev,[ :alpha => 0.78, :beta => 0.81] )	
		@fact ev.simMoments[:beta] => 0.81

		setValue(ev,4.2)
		@fact ev.value => 4.2
	end

	context("testing saving") do
		h5open("test5.h5", "w") do ff
			print("$ev")
			write(ff,"eval_test",ev)
			write(ff,"eval_list",[ev,ev2])
		end
		h5open("test5.h5", "r") do ff
			ev2 = readEval(ff,"eval_test")
			@fact ev2.status => ev.status
			evs = readEvalArray(ff,"eval_list")
			@fact evs[1].simMoments[:alpha] => ev.simMoments[:alpha]
			# readEval(ff,"eval_test",ev)
		end

	end
end



end # module 