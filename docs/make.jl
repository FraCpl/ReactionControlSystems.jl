using ReactionControlSystems
using Documenter

DocMeta.setdocmeta!(ReactionControlSystems, :DocTestSetup, :(using ReactionControlSystems); recursive=true)

makedocs(;
    modules=[ReactionControlSystems],
    authors="F. Capolupo",
    sitename="ReactionControlSystems.jl",
    format=Documenter.HTML(;
        canonical="https://FraCpl.github.io/ReactionControlSystems.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/FraCpl/ReactionControlSystems.jl",
    devbranch="master",
)
