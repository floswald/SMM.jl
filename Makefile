

all:
	julia -e 'include("src/MOpt.jl");include("./test/runtests.jl")'

testit:
	julia -e 'include("./test/runtests.jl")'


