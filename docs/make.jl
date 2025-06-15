using Documenter
using DocumenterVitepress

makedocs(;
    format=DocumenterVitepress.MarkdownVitepress(
        repo="github.com/YanisLeFur/Piqs.jl",
        devbranch="main",
        devurl="dev",
    ),
    sitename="Piqs.jl",
    authors="Yanis Le Fur",
    modules=[Piqs],
    pages=[
        "Home" => "index.md",
    ],
)


DocumenterVitepress.deploydocs(;
    repo="github.com/YanisLeFur/Piqs.jl",
    target=joinpath(@__DIR__, "build"),
    devbranch="main",
    branch="gh-pages",
    push_preview=true,
)
