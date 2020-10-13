using Documenter, ScenTrees

<<<<<<< HEAD
#const ASSETS = readdir(joinpath(@__DIR__, "src", "assets"))

makedocs(
	sitename =  "ScenTrees.jl",
	authors = "Kipngeno Kirui",
	clean = true,
	doctest = false,
	format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
	pages = ["Home" => "index.md",
		"Tutorials" => Any["tutorial/tutorial1.md",
				    "tutorial/tutorial2.md",
				    "tutorial/tutorial3.md",
				    "tutorial/tutorial31.md",
				    "tutorial/tutorial4.md",
				    "tutorial/tutorial41.md",
				    "tutorial/tutorial5.md"
				]
		]
)

deploydocs(
	repo = "github.com/kirui93/ScenTrees.jl.git",
	target = "build",
	deps = nothing,
	make = nothing
	#versions = ["stable" => "v^", "v#.#", "dev" => "master"]
=======
makedocs(
	sitename =  "ScenTrees.jl",
	authors = "Kipngeno Kirui",
	clean = false,
	doctest = true,
	format = Documenter.HTML(
		assets = ["exampleTree1.png"],
		prettyurls = get(ENV, "CI", nothing) == "true"),
	pages = [ "Home" => "index.md",
			"Tutorials" => "tutorial/tutorial1.md"],
)

deploydocs(deps = Deps.pip("mkdocs","python-markdown-math"),
	repo="github.com/kirui93/ScenTrees.jl.git"
>>>>>>> Deploy kirui93/ScenTrees.jl to github.com/kirui93/ScenTrees.jl.git:gh-pages
)
