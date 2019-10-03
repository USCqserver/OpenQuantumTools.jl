module QuantumAnnealingTools

using DocStringExtensions
using Reexport
using RecipesBase

import SpecialFunctions:trigamma
import Optim:optimize
import QuadGK:quadgk
import Arpack:eigs
import DiffEqBase:DEDataVector, DEDataMatrix, DEDataArray, ODEProblem, ODEFunction, DiscreteCallback, u_modified!, full_cache, solve, EnsembleSerial, EnsembleProblem, ODESolution

@reexport using LinearAlgebra
@reexport using SparseArrays
@reexport using QTBase
@reexport using LaTeXStrings

include("math_util.jl")

include("Bath/ohmic.jl")
include("Bath/hybridohmic.jl")
include("Bath/util.jl")

include("QSolver/util.jl")
include("QSolver/pausing_control.jl")
include("QSolver/schrodinger_solver.jl")
include("QSolver/unitary_solver.jl")
include("QSolver/von_neumann_solver.jl")
include("QSolver/ame_solver.jl")
include("QSolver/redfield_solver.jl")
include("QSolver/PTRE.jl")



include("plot_util/high_dpi.jl")
include("plot_util/hamiltonian_plot.jl")
include("plot_util/ode_sol.jl")
include("plot_util/bath_plot.jl")
include("plot_util/projected_system.jl")




export OhmicBath, Ohmic, γ, S, correlation, polaron_correlation, interpolate_spectral_density

export HybridOhmicBath, HybridOhmic, convolution_rate, Gₗ, Gₕ, half_width_half_maximum, bloch_rate, direct_integrate, spectrum_info, Sₕ, MRT_Γ

export info_freq, τ_SB, τ_B

export solve_unitary, solve_schrodinger, solve_von_neumann, solve_redfield, solve_ame, solve_af_rwa

export PausingControl, single_pausing

export SA_Δ², SA_redfield, SA_marcus, SA_Γ

export @publish, minimum_gap

end # end module