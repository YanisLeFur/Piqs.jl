using Documenter
using DocumenterVitepress
using PIQSolver

makedocs(;
    format=DocumenterVitepress.MarkdownVitepress(
        repo="github.com/YanisLeFur/PIQSolver.jl",
        devbranch="main",
        devurl="dev",
    ),
    sitename="PIQSolver.jl",
    authors="Yanis Le Fur",
    modules=[PIQSolver],
    checkdocs=:exports,
    pages=[
        "Home" => "index.md",
        "Getting started" => "getting_started.md",
        "Examples" => "examples.md",
        "API reference" => "api.md",
    ],
)


DocumenterVitepress.deploydocs(;
    repo="github.com/YanisLeFur/PIQSolver.jl",
    target=joinpath(@__DIR__, "build"),
    devbranch="main",
    branch="gh-pages",
    push_preview=true,
)
