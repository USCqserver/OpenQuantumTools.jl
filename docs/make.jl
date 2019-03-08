using Documenter, QTool

makedocs(sitename="QTool",
    pages = Any[
        "Home" => "index.md",
        "Library" => Any[
            "Public" => "lib/public.md"]
    ],
    format = Documenter.HTML(prettyurls = false)
)
