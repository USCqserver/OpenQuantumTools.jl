module QuantumAnnealingTools

using Reexport
import SpecialFunctions:trigamma
import Optim:optimize

@reexport using LinearAlgebra
@reexport using SparseArrays
@reexport using Arpack
@reexport using QInterpolations
@reexport using QIntegrations
@reexport using QTBase

include("QSolver/solvers.jl")
include("QSolver/me_util.jl")

include("Proj/proj_util.jl")
include("Proj/rotation_util.jl")

include("Bath/ohmic.jl")
include("Bath/hybridohmic.jl")

include("plot_util.jl")


export OhmicBath, Ohmic, γ, S, correlation, polaron_correlation, interpolate_spectral_density

export HybridOhmicBath, HybridOhmic, convolution_rate, GL, GH, half_width_half_maximum, tunneling_Δ, bloch_rate, spectrum_info

export  LowLevelSystem, RotatedTwoLevelSystem, proj_low_lvl, optimal_interaction_angle, get_dθ, rotate_lowlevel_system

export load_diff_eq
export calculate_unitary, check_unitary, solve_schrodinger, solve_von_neumann
export adiabatic_me_update!, solve_adiabatic_me

export @publish

end # end module
