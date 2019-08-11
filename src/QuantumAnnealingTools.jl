module QuantumAnnealingTools

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
include("Proj/proj_util.jl")
include("Proj/rotation_util.jl")

include("Bath/ohmic.jl")
include("Bath/hybridohmic.jl")

include("QSolver/util.jl")
include("QSolver/solvers.jl")
include("QSolver/schrodinger_solver.jl")

include("QSolver/ame_solvers.jl")


include("plot_util/high_dpi.jl")
include("plot_util/hamiltonian_plot.jl")
include("plot_util/ode_sol.jl")
include("plot_util/bath_plot.jl")




export OhmicBath, Ohmic, γ, S, correlation, polaron_correlation, interpolate_spectral_density

export HybridOhmicBath, HybridOhmic, convolution_rate, Gₗ, Gₕ, half_width_half_maximum, bloch_rate, direct_integrate, spectrum_info, Sₕ

export info_freq

export solve_unitary, solve_schrodinger, solve_von_neumann, solve_redfield, solve_davies

export  LowLevelSystem, RotatedTwoLevelSystem, proj_low_lvl, optimal_interaction_angle, get_dθ, rotate_lowlevel_system, @unitary_landau_zener, @unitary_interaction

export @publish, minimum_gap

end # end module
