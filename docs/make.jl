using Documenter, QuantumAnnealingTools

makedocs(sitename="QuantumAnnealingTools",
    pages = Any[
        "Home" => "index.md",
        "Tutorials" => Any[
            "Math Symbol" => "tutorials/math_symbol.md",
            "Hamiltonian" => "tutorials/hamiltonians.md",
            "Annealing" => "tutorials/annealing.md"
        ],
        "Library" => Any[
            "Base" => "lib/QTBase.md",
            #"Projection" => "lib/Proj.md",
            "Bath" => "lib/bath.md"
        ]
    ],
    format = Documenter.HTML(prettyurls = !("local" in ARGS))
)

deploydocs(
    repo = "github.com/USCqserver/QuantumAnnealingTools.jl.git",
)
