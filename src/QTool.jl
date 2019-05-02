module QTool

using Reexport
import SpecialFunctions:trigamma
import Optim:optimize

@reexport using LinearAlgebra
@reexport using SparseArrays
@reexport using Arpack

include("QInterpolate/QInterpolate.jl")
include("Integration/Integration.jl")
@reexport using .QInterpolate
@reexport using .Integration

include("QTBase/QTBase.jl")
@reexport using .QTBase
include("LinearOp/LinearOp.jl")
@reexport using .LinearOp
include("QSolver/QSolver.jl")
@reexport using.QSolver

include("Bath/ohmic.jl")
include("Bath/hybridohmic.jl")
include("Proj/proj_util.jl")

export OhmicBath, Ohmic, γ, S, correlation, polaron_correlation, interpolate_spectral_density

export HybridOhmicBath, HybridOhmic, convolution_rate, integral_rate, reorganization_energy, GL, GH, ohmic_correlation

export  LowLevelParams, RotatedTwoLevelParams, proj_low_lvl, optimal_interaction_angle, get_dθ, rotate_sys

end # end module
