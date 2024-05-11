using DeepFlex
using Documenter
ENV["GKSwstype"] = "100"
ENV["JULIA_DEBUG"] = "Documenter"

cp("./docs/Manifest.toml", "./docs/src/assets/Manifest.toml", force = true)
cp("./docs/Project.toml", "./docs/src/assets/Project.toml", force = true)

DocMeta.setdocmeta!(DeepFlex, :DocTestSetup,
    :(using DeepFlex); recursive = true)

makedocs(;
    modules = [DeepFlex],
    authors = "jing xin",
    sitename = "DeepFlex.jl",
    format = Documenter.HTML(assets = ["assets/favicon.ico"],
        canonical = "https://docs.sciml.ai/DeepFlex.jl/dev/"),
    clean = true,
    doctest = false,
    linkcheck = true,
    pages = [
        "Home" => "index.md",
        "Tutorials" => [
            "Run your first hydrological model" => "friction.md"
        ],
        "API" => "api.md"
    ]
)

deploydocs(;
    repo = "github.com/chooron/DeepFlex.jl",
    devbranch = "main",
    push_preview = true
)