using ClusterDepth
using Documenter
using Glob
using Literate
DocMeta.setdocmeta!(ClusterDepth, :DocTestSetup, :(using ClusterDepth); recursive=true)


GENERATED = joinpath(@__DIR__, "src")
for subfolder ∈ ["explanations","HowTo","tutorials","reference"]
    local SOURCE_FILES = Glob.glob(subfolder*"/*.jl", GENERATED)
    #config=Dict(:repo_root_path=>"https://github.com/unfoldtoolbox/UnfoldSim")
    foreach(fn -> Literate.markdown(fn, GENERATED*"/"*subfolder), SOURCE_FILES)

end
makedocs(;
    modules=[ClusterDepth],
    authors="Benedikt V. Ehinger",
    repo="https://github.com/s-ccs/ClusterDepth.jl/blob/{commit}{path}#{line}",
    sitename="ClusterDepth.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://s-ccs.github.io/ClusterDepth.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Demo" => "tutorials/demo.md"
    ],
)

deploydocs(;
    repo="github.com/s-ccs/ClusterDepth.jl",
    devbranch="main",
)
