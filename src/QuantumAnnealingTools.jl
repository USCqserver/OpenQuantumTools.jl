module QuantumAnnealingTools

using DocStringExtensions
using Reexport
using RecipesBase

import SpecialFunctions: trigamma
import Distributions: Exponential, product_distribution
import Optim: optimize
import QuadGK: quadgk
import Arpack: eigs
import DiffEqBase: DEDataVector,
                   DEDataMatrix,
                   DEDataArray,
                   ODEProblem,
                   ODEFunction,
                   DiscreteCallback,
                   ContinuousCallback,
                   u_modified!,
                   full_cache,
                   solve,
                   EnsembleSerial,
                   EnsembleProblem,
                   ODESolution,
                   DiffEqArrayOperator,
                   INITIALIZE_DEFAULT,
                   add_tstop!,
                   CallbackSet,
                   terminate!
import DiffEqCallbacks: FunctionCallingCallback
import DataStructures: compare
import NLSolversBase


@reexport using LinearAlgebra
@reexport using SparseArrays
@reexport using QTBase
@reexport using LaTeXStrings

include("math_util.jl")

include("Bath/custom_bath.jl")
include("Bath/ohmic.jl")
include("Bath/hybridohmic.jl")
include("Bath/onef.jl")
include("Bath/util.jl")

include("interaction.jl")

include("QControl/ame_trajectory_control.jl")
include("QControl/callback_lib.jl")
include("QControl/control_datatype.jl")
include("QControl/dd_control.jl")
include("QControl/onef_ode_control.jl")
include("QControl/manifold_projection.jl")

include("QSolver/util.jl")
include("QSolver/schrodinger_solver.jl")
include("QSolver/stochastic_schrodinger_solver.jl")
include("QSolver/unitary_solver.jl")
include("QSolver/von_neumann_solver.jl")
include("QSolver/ame_solver.jl")
include("QSolver/redfield_solver.jl")
include("QSolver/PTRE.jl")
include("QSolver/ensemble_builder.jl")

include("plot_util/high_dpi.jl")
include("plot_util/hamiltonian_plot.jl")
include("plot_util/ode_sol.jl")
include("plot_util/bath_plot.jl")
include("plot_util/projected_system.jl")


export OhmicBath,
       Ohmic,
       γ,
       S,
       correlation,
       polaron_correlation,
       interpolate_spectral_density,
       spectrum,
       CustomBath

export Interaction, InteractionSet

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

export EnsembleFluctuator

export info_freq, τ_SB, τ_B, coarse_grain_timescale

export solve_unitary,
       solve_schrodinger,
       solve_von_neumann,
       solve_redfield,
       solve_ame,
       solve_af_rwa,
       solve_stochastic_schrodinger,
       build_ensemble_problem,
       build_prob_func

export PausingControl, single_pausing, InstPulseControl, InstDEPulseControl, ControlSet, DEFAULT_INITIALIZER

export SA_Δ², SA_redfield, SA_marcus, SA_Γ, SA_τ, solve_SA, SA_lz_rotate

export @publish, minimum_gap, @highdpi

end # end module
