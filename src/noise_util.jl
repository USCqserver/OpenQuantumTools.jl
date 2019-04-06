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

function Base.show(io::IO, ::MIME"text/plain", m::OhmicBath)
    print(io, "Ohmic bath instance:\n", "η (unitless): ", m.η, "\n", "ωc (GHz): ", m.ωc/pi/2,
        "\n", "T (mK): ", beta_2_temperature(m.β))
end

"""
    HybridOhmic

A hybrid noise model with both low and high frequency noise. The high frequency noise is characterized by Ohmic bath and the low frequence noise is characterized by the MRT width `W`.

**Fields**
- `W` -- MRT width (2π GHz)
- `ϵ` -- reorganization energy ()
- `η` -- strength.
- `ωc` -- cutoff frequence.
- `β` -- inverse temperature.
"""
struct HybridOhmic
    W::Float64
    ϵ::Float64
    η::Float64
    ωc::Float64
    β::Float64
end

function Base.show(io::IO, ::MIME"text/plain", m::HybridOhmic)
    print(io, "Hybrid Ohmic bath instance:\n", "W (mK): ", freq_2_temperature(m.W/2/pi), "\n",
        "ϵ (GHz): ", m.ϵ/2/pi ,"\n",
        "η (unitless): ", m.η, "\n", "ωc (GHz): ", m.ωc/pi/2, "\n", "T (mK): ", beta_2_temperature(m.β))
end

function HybridOhmic(W, η, ωc, T)
    W = 2 * pi * temperature_2_freq(W)
    β = temperature_2_beta(T)
    ϵ = W^2 * β / 2
    ωc = 2 * π * ωc
    HybridOhmic(W, ϵ, η, ωc, β)
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
    correlation(τ, params::OhmicBath)

Calculate the correlation function of Ohmic bath.
"""
function correlation(τ, params::OhmicBath)
    x2 = 1 / params.β / params.ωc
    x1 = 1.0im * τ / params.β
    params.η * (trigamma(-x1+1+x2)+trigamma(x1+x2)) / params.β^2
end

function polaron_correlation(τ, params::Union{OhmicBath, HybridOhmic})
    nothing
end
