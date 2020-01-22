using Documenter, FinEtools, FinEtoolsDeforLinear, FinEtoolsVibInFluids

makedocs(
	modules = [FinEtoolsVibInFluids],
	doctest = false, clean = true,
	format = Documenter.HTML(prettyurls = false),
	authors = "Petr Krysl",
	sitename = "FinEtoolsVibInFluids.jl",
	pages = Any[
	"Home" => "index.md",
	"Guide" => "guide/guide.md",
	"Types and Functions" => Any[
		"man/types.md",
		"man/functions.md"]
		]
	)

deploydocs(
    repo = "github.com/PetrKryslUCSD/FinEtoolsVibInFluids.jl.git",
)
