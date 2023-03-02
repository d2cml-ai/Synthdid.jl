using Synthdid
using Documenter

DocMeta.setdocmeta!(Synthdid, :DocTestSetup, :(using Synthdid); recursive=true)

makedocs(;
    modules=[Synthdid],
    authors="tjhon <fr.jhonk@gmail.com> and contributors",
    repo="https://github.com/d2cml-ai/Synthdid.jl/blob/{commit}{path}#{line}",
    sitename="Synthdid.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://d2cml-ai.github.io/Synthdid.jl",
        edit_link="master",
        assets=String[]
    ),
    pages=[
        "Home" => "index.md",
    ]
)

deploydocs(;
    repo="github.com/d2cml-ai/Synthdid.jl",
    devbranch="master"
)
