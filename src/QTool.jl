module QTool

using Reexport
import SpecialFunctions:trigamma

@reexport using LinearAlgebra

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

export OhmicBath, Ohmic, γ, S, correlation, polaron_correlation, interpolate_spectral_density

export HybridOhmicBath, HybridOhmic, convolution_rate, integral_rate, reorganization_energy, GL, GH, ohmic_correlation

end # end module
