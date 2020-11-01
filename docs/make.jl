using Documenter, ScenTrees

#const ASSETS = readdir(joinpath(@__DIR__, "src", "assets"))

makedocs(
	sitename =  "ScenTrees.jl",
	authors = "Kipngeno Kirui",
<<<<<<< HEAD
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
)
=======
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

<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
if !Sys.iswindows()
	deploydocs(
		deps = Deps.pip("mkdocs","python-markdown-math"),
	   	repo = "github.com/kirui93/ScenTrees.jl.git",
	   	target = "site",
	   	make = () -> run(`mkdocs build`)  
	)
end
>>>>>>> e9b1bc9cdc5c989ee6e99a1505eeecf47d22e288
=======

deploydocs(repo = "github.com/kirui93/ScenTrees.jl.git")
>>>>>>> master
=======
if isCI
    deploydocs(repo = "github.com/kirui93/ScenTrees.jl.git")
end
>>>>>>> master
=======
#if isCI
=======
>>>>>>> master
deploydocs(
	repo = "github.com/kirui93/ScenTrees.jl.git",
	target = "build",
	versions = ["stable" => "v^", "v#.#", "dev" => "master"]
)
<<<<<<< HEAD
#end
>>>>>>> master
=======
>>>>>>> master
