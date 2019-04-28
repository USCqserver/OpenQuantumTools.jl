module QTool

using Reexport
import SpecialFunctions:trigamma

@reexport using LinearAlgebra

include("QTBase/QTBase.jl")
@reexport using .QTBase
include("LinearOp/LinearOp.jl")
@reexport using .LinearOp
include("QSolver/QSolver.jl")
@reexport using.QSolver

include("Bath/Bath.jl")

export OhmicBath, Ohmic, γ, S, correlation, polaron_correlation, interpolate_spectral_density

export HybridOhmicBath, HybridOhmic, convolution_rate, integral_rate, reorganization_energy, GL, GH, ohmic_correlation

end # end module
