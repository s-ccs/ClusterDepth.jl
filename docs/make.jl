using ClusterDepth
using Documenter

DocMeta.setdocmeta!(ClusterDepth, :DocTestSetup, :(using ClusterDepth); recursive=true)

makedocs(;
    modules=[ClusterDepth],
    authors="Benedikt V. Ehinger",
    repo="https://github.com/behinger/ClusterDepth.jl/blob/{commit}{path}#{line}",
    sitename="ClusterDepth.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://behinger.github.io/ClusterDepth.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/behinger/ClusterDepth.jl",
    devbranch="main",
)
