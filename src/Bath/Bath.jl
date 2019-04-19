module Bath

export OhmicBath, Ohmic, γ, S, correlation, polaron_correlation, interpolate_spectral_density
export HybridOhmicBath, HybridOhmic

include("../QTBase/QTBase.jl")
import .QTBase:beta_2_temperature, temperature_2_beta, freq_2_temperature, temperature_2_freq, cauchy_principal_value, construct_interpolations
import SpecialFunctions:trigamma
include("ohmic.jl")
include("hybridohmic.jl")

"""
    correlation(τ, params::OhmicBath)

Calculate the correlation function of Ohmic bath.
"""
function correlation(τ, params::Union{OhmicBath, HybridOhmicBath})
    x2 = 1 / params.β / params.ωc
    x1 = 1.0im * τ / params.β
    params.η * (trigamma(-x1+1+x2)+trigamma(x1+x2)) / params.β^2
end

end # end module
