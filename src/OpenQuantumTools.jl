module OpenQuantumTools

using DocStringExtensions
using Reexport
using RecipesBase

import SciMLBase:
    ODEProblem,
    EnsembleProblem,
    ODEFunction,
    DiscreteCallback,
    ContinuousCallback,
    u_modified!,
    full_cache,
    solve,
    ODESolution,
    DiffEqArrayOperator,
    CallbackSet,
    terminate!

import DiffEqCallbacks:
    FunctionCallingCallback, PresetTimeCallback, IterativeCallback

@reexport using LinearAlgebra
@reexport using SparseArrays
@reexport using OpenQuantumBase

include("math_util.jl")
include("QControl/callback_lib.jl")

include("QSolver/util.jl")
include("QSolver/closed_system_solvers.jl")
include("QSolver/lindblad_solver.jl")
include("QSolver/stochastic_schrodinger_solver.jl")
include("QSolver/ame_solver.jl")
include("QSolver/redfield_solver.jl")
include("QSolver/ensemble_builder.jl")

include("plot_util/hamiltonian_plot.jl")
include("plot_util/ode_sol.jl")
include("plot_util/bath_plot.jl")

export solve_unitary,
    solve_schrodinger,
    solve_von_neumann,
    solve_redfield,
    solve_ame,
    solve_cgme,
    solve_ule,
    solve_lindblad,
    build_ensembles

export DEFAULT_INITIALIZER
export InstPulseCallback, PositivityCheckCallback

end # end module
