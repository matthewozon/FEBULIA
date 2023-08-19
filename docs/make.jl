using Documenter, DocumenterTools
using FEBULIA

makedocs(
    sitename = "FEBULIA",
    format = Documenter.HTML(),
    modules = [FEBULIA],
    pages = [
    "Home" => "index.md",
    "PolyExp" => "PolyExp.md",
    "Inner products" => "inner_products.md",
    ],
    doctest = true,
)

deploydocs(repo = "github.com/matthewozon/FEBULIA.git",branch = "main") #, deploy_config=deploy_config_github) "github.com/matthewozon/FEBULIA.jl.git"




# pages = [
#     "Home" => "index.md",
#     "Manual" => Any[
#         "Guide" => "man/guide.md",
#         "man/syntax.md",
#     ],
#     ]

# deploy_config_github = Documenter.GitHubActions("github.com/matthewozon/FEBULIA.jl.git","push","refs/heads/main")