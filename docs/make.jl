using LXCat
using Documenter

DocMeta.setdocmeta!(LXCat, :DocTestSetup, :(using LXCat); recursive=true)

makedocs(;
    modules=[LXCat],
    authors="Jan Kuhfeld <jan.kuhfeld@rub.de> and contributors",
    repo="https://github.com/jqfeld/LXCat.jl/blob/{commit}{path}#{line}",
    sitename="LXCat.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://jqfeld.github.io/LXCat.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/jqfeld/LXCat.jl",
    devbranch="main",
)
