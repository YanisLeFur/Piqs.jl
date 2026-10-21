# Test groups, selected with the GROUP environment variable:
# - "Core" (default): unit tests and Aqua, run on every supported Julia version
# - "JET": static analysis with JET.jl, in its own environment (test/jet/Project.toml).
#   JET only supports specific Julia versions, so it runs in a separate CI job on a
#   pinned version instead of blocking the resolution of the main test environment.
const GROUP = get(ENV, "GROUP", "Core")

if GROUP == "Core"
    include("core.jl")
    include("aqua.jl")
elseif GROUP == "JET"
    using Pkg
    Pkg.activate(joinpath(@__DIR__, "jet"))
    Pkg.develop(PackageSpec(path=dirname(@__DIR__)))
    Pkg.instantiate()
    include("jet.jl")
else
    error("Unknown test GROUP = \"$GROUP\", expected \"Core\" or \"JET\".")
end
