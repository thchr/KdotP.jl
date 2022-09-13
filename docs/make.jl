using Documenter
using KdotP

# make sure we actually do `using KdotP` before calling doctests (as described in
# https://juliadocs.github.io/Documenter.jl/stable/man/doctests/#Module-level-metadata)
DocMeta.setdocmeta!(KdotP, :DocTestSetup, :(using KdotP); recursive=true)

makedocs(
    modules = [KdotP],
    sitename = "KdotP.jl",
    authors = "Thomas Christensen <tchr@mit.edu> and contributors",
    repo = "https://github.com/thchr/KdotP.jl/blob/{commit}{path}#L{line}",
    format=Documenter.HTML(;
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://thchr.github.io/KdotP.jl"
    ),
    pages = [
        "API" => "index.md"
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo   = "github.com/thchr/KdotP.jl.git",
    target = "build",
    deps   = nothing,
    make   = nothing,
    push_preview = true
)