using Documenter, QuantumAnnealingTools

makedocs(sitename="QuantumAnnealingTools",
    pages = Any[
        "Home" => "index.md",
        "Library" => Any[
            "Base" => "lib/QTBase.md",
            "Projection" => "lib/Proj.md",
            "Bath" => "lib/bath.md"
        ]
    ],
    format = Documenter.HTML(prettyurls = false)
)
