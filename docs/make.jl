using Documenter
using DocumenterVitepress
using Piqs

makedocs(;
    format=DocumenterVitepress.MarkdownVitepress(
        repo="github.com/YanisLeFur/Piqs.jl",
        devbranch="main",
        devurl="dev",
    ),
    sitename="Piqs.jl",
    authors="Yanis Le Fur",
    modules=[Piqs],
    checkdocs=:exports,
    pages=[
        "Home" => "index.md",
        "Getting started" => "getting_started.md",
        "Examples" => "examples.md",
        "API reference" => "api.md",
    ],
)


DocumenterVitepress.deploydocs(;
    repo="github.com/YanisLeFur/Piqs.jl",
    target=joinpath(@__DIR__, "build"),
    devbranch="main",
    branch="gh-pages",
    push_preview=true,
)
