using Piqs
using Documenter
using DocumenterVitepress

makedocs(;
    modules=[Piqs],
    authors="Yanis Le Fur",
    sitename="Piqs.jl",
    repo=Remotes.GitHub("https://github.com/YanisLeFur/Piqs.jl"),
    format=DocumenterVitepress.MarkdownVitepress(
        repo="github.com/YanisLeFur/Piqs.jl",
        devbranch="main",
        devurl="dev",
    ),
    pages=[
        "Home" => "index.md",
    ],
    draft=false,
    doctest=true,
    plugins=[],
)

DocumenterVitepress.deploydocs(;
    repo="github.com/YanisLeFur/Piqs.jl",
    target=joinpath(@__DIR__, "build"),
    devbranch="main",
    branch="gh-pages",
    push_preview=true,
)
