module QSolver

import LinearAlgebra:mul!, Matrix, I, Diagonal, axpy!
include("math_util.jl")
include("me_util.jl")
include("solvers.jl")

export load_diff_eq
export calculate_unitary, check_unitary, solve_schrodinger, solve_von_neumann
export adiabatic_me_update!, solve_adiabatic_me

function load_diff_eq()
    @eval QSolver begin
        using DifferentialEquations:ODEProblem, solve, Tsit5, ODEFunction
    end
end

end  # module QSolver
