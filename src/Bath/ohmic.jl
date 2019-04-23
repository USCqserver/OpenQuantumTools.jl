"""
    OhmicBath

Ohmic bath object to hold a particular parameter set.

**Fields**
- `η` -- strength.
- `ωc` -- cutoff frequence.
- `β` -- inverse temperature.
"""
struct OhmicBath
    η::Float64
    ωc::Float64
    β::Float64
end

"""
    Ohmic(η, fc, T)

Construct OhmicBath from parameters with physical unit: `η`--unitless interaction strength; `fc`--cutoff frequency in GHz; `T`--temperature in mK.
"""
function Ohmic(η, fc, T)
    ωc = 2 * π * fc
    β = temperature_2_beta(T)
    OhmicBath(η, ωc, β)
end

function Base.show(io::IO, ::MIME"text/plain", m::OhmicBath)
    print(io, "Ohmic bath instance:\n", "η (unitless): ", m.η, "\n", "ωc (GHz): ", m.ωc/pi/2,
        "\n", "T (mK): ", beta_2_temperature(m.β))
end

"""
    γ(w::Float64, params::OhmicBath)

Calculate real part of Ohmic spectrum. The value between [-1000, 1000]*eps() is set to ``2πη/β``.
"""
function γ(w::Float64, params::OhmicBath)
    if w > 1000 * eps()
        return 2 * pi * params.η * w * exp(-w/params.ωc) / (1 - exp(-params.β*w))
    elseif w < - 1000 * eps()
        temp = exp(params.β * w)
        return -2 * pi * params.η * w * exp(w/params.ωc) * temp / (1 - temp)
    else
        return 2* pi* params.η / params.β
    end
end

"""
    S(w::Float64, params::OhmicBath; atol=1e-7)

Calculate the Lamb shift of Ohmic spectrum. `atol` is the absolute tolerance for Cauchy principal value integral.
"""
function S(w::Float64, params::OhmicBath; atol=1e-7)
    f(x)= γ(x, params)
    res = cauchy_principal_value(f, w, atol=atol)
    -res[1]/2/pi
end

"""
    polaron_correlation(τ, params::OhmicBath)

Calculate the polaron transformed correlation function of Ohmic bath.
"""
function polaron_correlation(τ, params::OhmicBath)
    res = (1+1.0im*params.ωc*τ)^(-4*params.η)
    if !isapprox(τ, 0, atol = 1e-9)
        x = π * τ / params.β
        res *= ( x / sinh(x) )^(4 * params.η)
    end
    res
end

"""
    interpolate_spectral_density(ω_grid, params::OhmicBath)

Calculate the Ohmic bath spectral density S on grid `ω_grid`, and construct interpolation objects for it. A separate function for γ is also returned without doing interpolation.
"""
function interpolate_spectral_density(ω_grid::AbstractRange{T}, params::OhmicBath) where T<: Number
    s_list = [S(ω, params) for ω in ω_grid]
    s_itp = construct_interpolations(ω_grid, s_list)
    (ω)->γ(ω, params), s_itp
end
