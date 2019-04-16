"""
    HybridOhmicBath

A hybrid noise model with both low and high frequency noise. The high frequency noise is characterized by Ohmic bath and the low frequence noise is characterized by the MRT width `W`.

**Fields**
- `W` -- MRT width (2π GHz)
- `ϵ` -- reorganization energy (2π GHz)
- `η` -- strength.
- `ωc` -- cutoff frequence.
- `β` -- inverse temperature.
"""
struct HybridOhmicBath
    W::Float64
    ϵ::Float64
    η::Float64
    ωc::Float64
    β::Float64
end

function polaron_correlation(τ, params::HybridOhmicBath)
    ohmic_part = (1+1.0im*params.ωc*τ)^(-4*params.η)
    if !isapprox(τ, 0, atol = 1e-9)
        x = π * τ / params.β
        ohmic_part *= ( x / sinh(x) )^(4 * params.η)
    end
    slow_part = exp( - 2.0 * params.W^2 * τ^2 - 4.0im * τ * params.ϵ)
    ohmic_part * slow_part
end
