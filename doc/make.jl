using Documenter
using HydroModels

makedocs(
    sitename = "HydroModels.jl",
    authors = "xin jing",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://chooron.github.io/HydroModels.jl",
        assets = ["assets/icons.ico"],
        highlights = ["yaml"],
        collapselevel = 2,
        sidebar_sitename = false
    ),
    modules = [HydroModels],
    clean = true,
    doctest = true,
    checkdocs = :none,
    pages = [
        "Home" => "index.md",
        # "Getting Started" => "getting_started.md",
        "Tutorials" => [
            "Run a Bucket Model" => "tutorials/run_a_bucket.md",
            "Run ExpHydro Model" => "tutorials/run_a_exphydro_model.md",
            # "Spatial Model Example" => "tutorials/spatial_model.md"
        ],
        # "API Reference" => [
        #     "Model Types" => "api/types.md",
        #     "Model Functions" => "api/functions.md",
        #     "Utility Functions" => "api/utils.md"
        # ],
        # "Contributing" => "contributing.md"
    ]
)

deploydocs(
    repo = "github.com/chooron/HydroModels.jl",
    devbranch = "main",
    push_preview = true
)