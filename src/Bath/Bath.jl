module Bath

export OhmicBath, Ohmic, γ, S, correlation, polaron_correlation, interpolate_spectral_density
export HybridOhmicBath, HybridOhmic

include("../QTBase/QTBase.jl")
import .QTBase:beta_2_temperature, temperature_2_beta, freq_2_temperature, temperature_2_freq, cauchy_principal_value, construct_interpolations, RotatedTwoLevelParams
import SpecialFunctions:trigamma
include("ohmic.jl")
include("hybridohmic.jl")

end # end module
