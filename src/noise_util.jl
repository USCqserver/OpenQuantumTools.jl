"""
    OhmicBath

Ohmic bath object to hold a particular parameter set.

**Fields**
- `η` -- strength.
- `wc` -- cutoff frequence.
- `β` -- inverse temperature.
"""
struct OhmicBath
    η::Float64
    wc::Float64
    β::Float64
end

"""
    γ(w::Float64, params::OhmicBath)

Calculate real part of Ohmic spectrum. The value between [-1000, 1000]*eps() is set to ``2πη/β``.
"""
function γ(w::Float64, params::OhmicBath)
    if w > 1000 * eps()
        return 2 * pi * params.η * w * exp(-w/params.wc) / (1 - exp(-params.β*w))
    elseif w < - 1000 * eps()
        temp = exp(params.β * w)
        return -2 * pi * params.η * w * exp(w/params.wc) * temp / (1 - temp)
    else
        return 2* pi* params.η / params.β
    end
end

"""
    correlation(τ, params::OhmicBath)

Calculate the correlation function of Ohmic bath.
"""
function correlation(τ, params::OhmicBath)
    x2 = 1 / params.β / params.wc
    x1 = 1.0im * τ / params.β
    params.η * (trigamma(-x1+1+x2)+trigamma(x1+x2)) / params.β^2
end
