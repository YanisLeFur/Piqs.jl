using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using Documenter
using DocumenterVitepress
using Piqs

const PAGES = ["Home" => "index.md"]

makedocs(;
    modules=[
        QuantumToolbox,
        Base.get_extension(QuantumToolbox, :QuantumToolboxMakieExt),
    ],
    authors="Yanis Le Fur",
    repo=Remotes.GitHub("Piqs.jl"),
    sitename="Piqs.jl",
    pages=PAGES,
    format=DocumenterVitepress.MarkdownVitepress(
        repo="github.com/YanisLeFur/Piqs.jl",
        devbranch="main",
        devurl="dev",
    ),
)


DocumenterVitepress.deploydocs(;
    repo="github.com/YanisLeFur/Piqs.jl",
    target=joinpath(@__DIR__, "build"),
    devbranch="main",
    branch="gh-pages",
    push_preview=true,
)