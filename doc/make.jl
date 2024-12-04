push!(LOAD_PATH, "../src/")
using Documenter
using HydroModels

makedocs(
    sitename = "HydroModels.jl",
    authors = "xin jing",
    format = Documenter.HTML(
        # 启用 pretty URLs，移除 .html 后缀
        prettyurls = get(ENV, "CI", nothing) == "true",
        # 设置文档的规范 URL
        canonical = "https://chooron.github.io/HydroModels.jl",
        # 设置资源文件
        assets = ["assets/icons.ico"],
        # 配置侧边栏
        collapselevel = 2,
        sidebar_sitename = true,
        # # 添加 favicon
        # favicon = "assets/icons.ico"
    ),
    # 配置模块
    modules = [HydroModels],
    # 配置页面结构
    pages = [
        "Home" => "index.md",
        "Tutorials" => [
            "Run a Bucket Model" => "tutorials/run_a_bucket.md",
            "Run ExpHydro Model" => "tutorials/run_a_exphydro_model.md"
        ]
    ],
    # 其他选项
    checkdocs = :none,
    clean = true
)

# 部署配置
deploydocs(
    repo = "github.com/chooron/HydroModels.jl",
    devbranch = "main",
    push_preview = true
)