

# TESTING Eval 
# ============




@testset "Testing Eval object" begin

	p    = Dict(:a => 3.1 , :b => 4.9)
	moms = DataFrame(name=[:alpha;:beta;:gamma],value=[0.8;0.7;0.5],weight=rand(3))

	p2    = Dict(:a => 3.1 , :b => 4.9)
	moms2 = DataFrame(name=[:alpha;:beta;:gamma],value=[0.8;0.7;0.5],weight=[0.1;0.4;0.5])

	p3    = MomentOpt.OrderedDict(:a => 3.1 , :b => 4.9)
	moms3 = DataFrame(name=[:alpha;:beta;:gamma],value=[0.8;0.7;0.5],weight=[0.1;0.4;0.5])

	ev = Eval(p,moms)
	setMoment(ev,moms)
	ev.status=1
	ev2= Eval(p2,moms2)	
	setMoment(ev2,moms2)
	ev2.status=1
	ev3= Eval(p3,moms3)	
	# setMoment(ev3,moms3)
	ev3.status=1

	@testset "testing param" begin

		@test param(ev,:a) == 3.1
		@test param(ev,[:a,:b]) == [3.1;4.9]

		@test param(ev2,:a) == 3.1
		@test param(ev2,[:a;:b]) == [3.1;4.9]
		@test paramd(ev2) == Dict(:a => 3.1 , :b => 4.9)

		@test param(ev3,:a) == 3.1
		@test param(ev3,[:a;:b]) == [3.1;4.9]
		@test paramd(ev3) == Dict(:a => 3.1 , :b => 4.9)

		myp = MyP()
		MomentOpt.fill(myp,ev)
		@test myp.a == 3.1
	end

	@testset "testing moments" begin
		
		@test dataMoment(ev2,:alpha) == 0.8
		@test dataMomentW(ev2,:alpha) == 0.1	

		setMoment(ev,:alpha,0.78)	
		setMoment(ev2,DataFrame(name="alpha" , value = 0.78))	

		@test ev.simMoments[:alpha] == 0.78
		@test ev2.simMoments[:alpha] == 0.78

		setMoment(ev,Dict( :alpha => 0.78, :beta => 0.81) )	
		@test ev.simMoments[:beta] == 0.81

		MomentOpt.setValue(ev,4.2)
		@test ev.value == 4.2
	end

	# @testset "testing saving/loading" begin
	# 	MomentOpt.FileIO.save("test5.jld2",Dict("ev"=>ev,"evs"=>[ev,ev2]))
	# 	sleep(1)
	# 	ev2_ = MomentOpt.FileIO.load("test5.jld2", "ev")
	# 	@test ev2_ == ev

	# 	evs = MomentOpt.FileIO.load("test5.jld2","evs")
	# 	@test evs[1] == ev
	# 	@test evs[2] == ev2
	# 	rm("test5.jld2")

	# end
end



