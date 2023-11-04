#!/usr/bin/env -S julia --startup-file=no

using Documenter
using KangarooTwelve
using Org

orgfiles = filter(f -> endswith(f, ".org"),
                  readdir(joinpath(@__DIR__, "src"), join=true))

for orgfile in orgfiles
    mdfile = replace(orgfile, r"\.org$" => ".md")
    read(orgfile, String) |>
        c -> Org.parse(OrgDoc, c) |>
        o -> sprint(markdown, o) |>
        s -> replace(s, r"\.org]" => ".md]") |>
        m -> write(mdfile, m)
end

makedocs(;
    modules=[KangarooTwelve],
    format=Documenter.HTML(),
    pages=[
        "KangarooTwelve" => "index.md",
        "Internals" => "internals.md",
    ],
    repo="https://github.com/tecosaur/KangarooTwelve.jl/blob/{commit}{path}#L{line}",
    sitename="KangarooTwelve.jl",
    authors = "tecosaur and contributors: https://github.com/tecosaur/KangarooTwelve.jl/graphs/contributors",
    warnonly = [:missing_docs],
)

deploydocs(repo="github.com/tecosaur/KangarooTwelve.jl")
