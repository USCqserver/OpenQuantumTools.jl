module Bath

export OhmicBath, γ, S, correlation, polaron_correlation
export HybridOhmicBath

include("../Integration/Integration.jl")
import .Integration:cauchy_principal_value
import SpecialFunctions:trigamma
include("ohmic.jl")
include("hybridohmic.jl")

end # end module
