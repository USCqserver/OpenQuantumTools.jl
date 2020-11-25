using Documenter, OpenQuantumTools
using DocumenterTools: Themes

Themes.compile(joinpath(@__DIR__, "src/assets/design.scss"))

makedocs(
    format=Documenter.HTML(
        prettyurls=false,
        assets = ["assets/design.css"],
        mathengine=MathJax(),
    ),
    sitename="OpenQuantumTools",
    pages = Any[
        "Home" => "index.md",
        "Tutorials" => Any[
            "Math Symbol" => "tutorials/math_symbol.md",
            "Hamiltonian" => "tutorials/hamiltonians.md",
            "Coupling" => "tutorials/couplings.md",
            "Bath" => "tutorials/bath.md",
            "Annealing" => "tutorials/annealing.md",
            "Solver" => "tutorials/solver.md"
        ],
        "Library" => Any[
            "Hamiltonian" => "lib/hamiltonian.md",
            "Matrix Utility" => "lib/matrix_util.md",
            "Coupling" => "lib/coupling.md",
            "Bath"=>"lib/bath.md",
            "Interaction"=>"lib/interactions.md",
            "Annealing"=>"lib/annealing.md",
            "Solver"=>"lib/solvers.md",
            "Callback"=>"lib/callbacks.md"
        ]
    ]
)

deploydocs(
    repo = "github.com/USCqserver/OpenQuantumTools.jl.git",
)
