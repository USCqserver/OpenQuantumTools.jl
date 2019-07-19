module QuantumAnnealingTools

using Reexport
import SpecialFunctions:trigamma
import Optim:optimize
import QuadGK:quadgk
import Arpack:eigs
import DiffEqBase:DEDataVector, DEDataMatrix, DEDataArray, ODEProblem, ODEFunction, DiscreteCallback, u_modified!, full_cache, solve

@reexport using LinearAlgebra
@reexport using SparseArrays
@reexport using QTBase

include("Proj/proj_util.jl")
include("Proj/rotation_util.jl")

include("Bath/ohmic.jl")
include("Bath/hybridohmic.jl")

include("QSolver/solvers.jl")


include("plot_util.jl")


export OhmicBath, Ohmic, γ, S, correlation, polaron_correlation, interpolate_spectral_density

export HybridOhmicBath, HybridOhmic, convolution_rate, Gₗ, Gₕ, half_width_half_maximum, bloch_rate, direct_integrate, spectrum_info, Sₕ

export solve_unitary, solve_schrodinger, solve_von_neumann, solve_redfield, solve_davies

export  LowLevelSystem, RotatedTwoLevelSystem, proj_low_lvl, optimal_interaction_angle, get_dθ, rotate_lowlevel_system, @unitary_landau_zener, @unitary_interaction

export @publish

end # end module
