using Documenter, QTool

makedocs(sitename="QTool",
    pages = Any[
        "Home" => "index.md",
        "Library" => Any[
            "Base" => "lib/QTBase.md",
            "Bath" => "lib/bath.md"
        ]
    ],
    format = Documenter.HTML(prettyurls = false)
)
