module QuantumAnnealingTools

using Reexport
import SpecialFunctions:trigamma
import Optim:optimize
import QuadGK:quadgk

@reexport using LinearAlgebra
@reexport using SparseArrays
@reexport using Arpack
@reexport using QInterpolations
@reexport using QTBase

include("QSolver/solvers.jl")
include("QSolver/me_util.jl")

include("Proj/proj_util.jl")
include("Proj/rotation_util.jl")

include("Bath/ohmic.jl")
include("Bath/hybridohmic.jl")

include("plot_util.jl")


export OhmicBath, Ohmic, γ, S, correlation, polaron_correlation, interpolate_spectral_density

export HybridOhmicBath, HybridOhmic, convolution_rate, Gₗ, Gₕ, half_width_half_maximum, bloch_rate, direct_integrate, spectrum_info, Sₕ

export  LowLevelSystem, RotatedTwoLevelSystem, proj_low_lvl, optimal_interaction_angle, get_dθ, rotate_lowlevel_system, @unitary_landau_zener, @unitary_interaction

export load_diff_eq
export calculate_unitary, check_unitary, solve_schrodinger, solve_von_neumann
export adiabatic_me_update!, solve_adiabatic_me

export @publish

end # end module
