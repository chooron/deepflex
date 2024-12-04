using Documenter
using HydroModels

makedocs(
    sitename="HydroModels.jl",
    authors="xin jing",
    format=Documenter.HTML(assets=["assets/icons.ico"]),
    modules=[HydroModels],
    checkdocs=:none,
    pages=[
        "Home" => "index.md",
        "Tutorials" => [
            "Run a Bucket Model" => "tutorials/run_a_bucket.md",
            "Run ExpHydro Model" => "tutorials/run_a_exphydro_model.md"
        ]
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=
deploydocs(
    repo = "github.com/USERNAME/HydroModels.jl.git"
)
=#