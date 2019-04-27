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

function polaron_correlation(τ, a, params::HybridOhmicBath)
    η = a * params.η
    ϵ = a * params.ϵ
    W2 = a * params.W^2
    ohmic_part = (1+1.0im*params.ωc*τ)^(-4*η)
    if !isapprox(τ, 0, atol = 1e-9)
        x = π * τ / params.β
        ohmic_part *= ( x / sinh(x) )^(4 * η)
    end
    slow_part = exp( - 2.0 * W2 * τ^2 - 4.0im * τ * ϵ)
    ohmic_part * slow_part
end

function HybridOhmic(W, η, fc, T)
    W = 2 * pi * temperature_2_freq(W)
    β = temperature_2_beta(T)
    ϵ = W^2 * β / 2
    ωc = 2 * π * fc
    HybridOhmicBath(W, ϵ, η, ωc, β)
end

function Base.show(io::IO, ::MIME"text/plain", m::HybridOhmicBath)
    print(io, "Hybrid Ohmic bath instance:\n", "W (mK): ", freq_2_temperature(m.W/2/pi), "\n",
        "ϵ (GHz): ", m.ϵ/2/pi ,"\n",
        "η (unitless): ", m.η, "\n", "ωc (GHz): ", m.ωc/pi/2, "\n", "T (mK): ", beta_2_temperature(m.β))
end

function γ(w::Float64, params::HybridOhmicBath)
    if w > 1000 * eps()
        return 8 * pi * params.η * w * exp(-w/params.ωc) / (1 - exp(-params.β*w))
    elseif w < - 1000 * eps()
        temp = exp(params.β * w)
        return -8 * pi * params.η * w * exp(w/params.ωc) * temp / (1 - temp)
    else
        return 8* pi* params.η / params.β
    end
end

function correlation(τ, params::HybridOhmicBath)
    x2 = 1 / params.β / params.ωc
    x1 = 1.0im * τ / params.β
    params.η * (trigamma(-x1+1+x2)+trigamma(x1+x2)) / params.β^2
end

function convolution_rate(sys, bath::HybridOhmicBath)
    Γ10 = []
    Γ01 = []
    for i in eachindex(sys.s)
        T_bar = sys.T[i] - sys.d[i] * bath.ϵ
        A = abs2(T_bar - sys.ω[i] * sys.c[i] / sys.a[i])
        B = (sys.a[i] * sys.b[i] - abs2(sys.c[i])) / sys.a[i]^2
        W = sys.a[i] * bath.W^2
        ϵ = sys.a[i] * bath.ϵ
        γ2 = (sys.a[i] * γ(0.0, bath) / 2)^2
        Δ = (ω) -> A + B * (ω^2 + W)
        GL = (ω) -> sqrt(2*π/W) * exp(-(ω-ϵ)^2/2/W)
        GH = (ω) -> sys.a[i] * γ(ω, bath) / (ω^2 + γ2)
        integrand_12 = (ω)->Δ(ω) * GL(sys.ω[i] - ω) * GH(ω)
        integrand_21 = (ω)->Δ(ω) * GL(-sys.ω[i] - ω) * GH(ω)
        push!(Γ10, integrate_1d(integrand_12, -Inf, Inf)[1])
        push!(Γ01, integrate_1d(integrand_21, -Inf, Inf)[1])
    end
    Γ10/2/pi, Γ01/2/pi
end

function integral_rate(sys, bath::HybridOhmicBath)
    Γ10 = []
    Γ01 = []
    for i in eachindex(sys.s)
        T_bar = sys.T[i] - (sys.d[i] + sys.c[i]) * bath.ϵ
        integrand_12 = (x)->(sys.b[i] * correlation(x, bath) + abs2(T_bar)) * polaron_correlation(x, sys.a[i], bath) * exp(1.0im * sys.ω[i] * x)
        integrand_21 = (x)->(sys.b[i] * correlation(x, bath) + abs2(T_bar)) * polaron_correlation(x, sys.a[i], bath) * exp(-1.0im * sys.ω[i] * x)
        push!(Γ10, integrate_1d(integrand_12, -Inf, Inf)[1])
        push!(Γ01, integrate_1d(integrand_21, -Inf, Inf)[1])
    end
    Γ10, Γ01
end
