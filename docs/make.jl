using Documenter, QuantumAnnealingTools

makedocs(
    format=Documenter.HTML(
        prettyurls=false,
        mathengine=MathJax(),
    ),
    sitename="QuantumAnnealingTools",
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
            "Matrix Utilities" => "lib/matrix_util.md",
            "Coupling" => "lib/coupling.md",
            "Bath"=>"lib/bath.md",
            "Annealing"=>"lib/annealing.md",
            "Solvers"=>"lib/solvers.md"
        ]
    ]
)

deploydocs(
    repo = "github.com/USCqserver/QuantumAnnealingTools.jl.git",
)
