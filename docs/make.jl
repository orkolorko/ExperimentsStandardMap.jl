using ExperimentsStandardMap
using Documenter

DocMeta.setdocmeta!(ExperimentsStandardMap, :DocTestSetup, :(using ExperimentsStandardMap); recursive=true)

makedocs(;
    modules=[ExperimentsStandardMap],
    authors="Isaia Nisoli nisoli@im.ufrj.br and contributors",
    repo="https://github.com/orkolorko/ExperimentsStandardMap.jl/blob/{commit}{path}#{line}",
    sitename="ExperimentsStandardMap.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://orkolorko.github.io/ExperimentsStandardMap.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/orkolorko/ExperimentsStandardMap.jl",
    devbranch="main",
)
