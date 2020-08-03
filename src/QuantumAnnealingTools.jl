module QuantumAnnealingTools

using DocStringExtensions
using Reexport
using RecipesBase

import Optim: optimize
import DiffEqBase:
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
@reexport using QTBase

include("math_util.jl")
include("QControl/callback_lib.jl")

#include("QControl/ame_trajectory_control.jl")
#include("QControl/control_datatype.jl")
include("QSolver/util.jl")
include("QSolver/interaction_dispatch.jl")
include("QSolver/closed_system_solvers.jl")
include("QSolver/stochastic_schrodinger_solver.jl")
include("QSolver/ame_solver.jl")
include("QSolver/redfield_solver.jl")
#include("QSolver/PTRE.jl")
include("QSolver/ensemble_builder.jl")

include("plot_util/hamiltonian_plot.jl")
include("plot_util/ode_sol.jl")
include("plot_util/bath_plot.jl")
include("plot_util/projected_system.jl")

export HybridOhmicBath,
    HybridOhmic,
    convolution_rate,
    Gₗ,
    Gₕ,
    half_width_half_maximum,
    bloch_rate,
    direct_integrate,
    spectrum_info,
    Sₕ,
    MRT_Γ

export solve_unitary,
    solve_schrodinger,
    solve_von_neumann,
    solve_redfield,
    solve_ame,
    solve_CGME,
    build_ensembles

export DEFAULT_INITIALIZER
export minimum_gap
export InstPulseCallback, PositivityCheckCallback


#export SA_Δ², SA_redfield, SA_marcus, SA_Γ, SA_τ, solve_SA, SA_lz_rotate

end # end module
