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
    γ(ω, params::OhmicBath)

Calculate Ohmic spectrum density, defined as a full Fourier transform on the bath correlation function.
"""
function γ(ω, params::OhmicBath)
    if isapprox(ω, 0.0, atol = 1e-9)
        return 2* pi* params.η / params.β
    else
        return 2 * pi * params.η * ω * exp(-abs(ω)/params.ωc) / (1 - exp(-params.β*ω))
    end
end

"""
    S(w, params::OhmicBath; atol=1e-7)

Calculate the Lamb shift of Ohmic spectrum. `atol` is the absolute tolerance for Cauchy principal value integral.
"""
function S(w, params::OhmicBath; atol=1e-7)
    f(x)= γ(x, params)
    g(x) = f(x) / (x - w)
    cpv, cperr = cpvagk(f, w, w-1.0, w+1.0)
    negv, negerr = quadgk(g, -Inf, w-1.0)
    posv, poserr = quadgk(g, w+1.0, Inf)
    v = cpv + negv + posv
    err = cperr + negerr + poserr
    if (err > atol) || (isnan(err))
        @warn "Absolute error of integration is larger than the tolerance."
    end
    -v/2/pi
end

"""
    correlation(τ, params::OhmicBath)

Calculate the correlation function of Ohmic bath.
"""
function correlation(τ, params::OhmicBath)
    x2 = 1 / params.β / params.ωc
    x1 = 1.0im * τ / params.β
    params.η * (trigamma(-x1+1+x2)+trigamma(x1+x2)) / params.β^2
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
    interpolate_spectral_density(ω_grid::AbstractRange{T}, params::OhmicBath) where T<:Number

Calculate the Ohmic bath spectral density S on grid `ω_grid`, and construct interpolation objects for it. A separate function for γ is also returned without doing interpolation.
"""
function interpolate_spectral_density(ω_grid::AbstractRange{T}, params::OhmicBath) where T<:Number
    s_list = [S(ω, params) for ω in ω_grid]
    s_itp = construct_interpolations(ω_grid, s_list)
    (ω)->γ(ω, params), s_itp
end
