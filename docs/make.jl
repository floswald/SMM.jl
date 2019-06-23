
using Documenter
using MomentOpt

makedocs(
    format = Documenter.HTML(
        prettyurls = false,
        canonical = "https://floswald.github.io/MomentOpt.jl/stable",
        analytics = "UA-41584331-5"),
    sitename = "MomentOpt.jl",
    modules = [MomentOpt],
    pages = [
    	"Introduction" => "index.md",
    	"Eval Object" => "eval.md",
    	"Slices" => "slices.md",
    	"Moment Algorithms" => "algo.md",
        "Econometrics" => "metrics.md",
    	"Examples" => "examples.md"]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/floswald/MomentOpt.jl.git"
)
