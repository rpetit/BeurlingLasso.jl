using Documenter, BeurlingLasso

makedocs(
    modules = [BeurlingLasso],
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", nothing) == "true"),
    authors = "Romain Petit",
    sitename = "BeurlingLasso",
    pages = [
        "Introduction" => "index.md",
        "Examples" => "examples.md"
    ]
    # strict = true,
    # clean = true,
    # checkdocs = :exports,
)

deploydocs(
    repo = "github.com/rpetit/BeurlingLasso.jl.git",
)
