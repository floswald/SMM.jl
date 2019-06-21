
using Documenter
using MomentOpt

makedocs(
    format = Documenter.HTML(),
    sitename = "MomentOpt.jl",
    modules = [MomentOpt],
    pages = [
    	"Getting Started" => "index.md",
    	"Eval Object" => "eval.md",
    	"Slices" => "slices.md",
    	"Moment Algorithms" => "algo.md",
    	"Examples" => "examples.md"]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/floswald/MomentOpt.jl.git",
    devbranch = "v0.7.2"
)
