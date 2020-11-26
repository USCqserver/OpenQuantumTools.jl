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
        "Basics" => Any[
            "Math Symbol" => "basics/math_symbol.md",
            "Hamiltonian" => "basics/hamiltonians.md",
            "Coupling" => "basics/couplings.md",
            "Bath" => "basics/bath.md",
            "Annealing" => "basics/annealing.md",
            "Solver" => "basics/solver.md"
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
