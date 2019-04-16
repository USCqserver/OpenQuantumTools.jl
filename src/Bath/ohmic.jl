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
