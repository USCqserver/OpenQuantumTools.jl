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

"""
    HybridOhmic(W, η, fc, T)

Construct HybridOhmicBath object with parameters in physical units. W: MRT width (mK); η: interaction strength (unitless); fc: Ohmic cutoff frequency (GHz); T: temperature (mK).
"""
function HybridOhmic(W, η, fc, T)
    W = 2 * pi * temperature_2_freq(W)
    β = temperature_2_beta(T)
    ϵ = W^2 * β / 2
    ωc = 2 * π * fc
    HybridOhmicBath(W, ϵ, η, ωc, β)
end

function correlation(τ, bath::HybridOhmicBath, a=1)
    W² = 4 * a * bath.W^2
    x2 = 1 / bath.β / bath.ωc
    x1 = 1.0im * τ / bath.β
    W² + 4 * bath.η * (trigamma(-x1+1+x2)+trigamma(x1+x2)) / bath.β^2
end

"""
    polaron_correlation(τ, bath::HybridOhmicBath[, a=1])

Calculate polaron transformed correlation function of HybridOhmicBath 'bath' at time 'τ' with relative strength `a`. The effective strength will be `` a * W^2`` and `` a * η ``.
"""
function polaron_correlation(τ, bath::HybridOhmicBath, a=1)
    η = a * bath.η
    ϵ = a * bath.ϵ
    W² = a * bath.W^2
    ohmic_part = (1+1.0im*bath.ωc*τ)^(-4*η)
    if !isapprox(τ, 0, atol = 1e-9)
        x = π * τ / bath.β
        ohmic_part *= ( x / sinh(x) )^(4*η)
    end
    slow_part = exp( - 2.0 * W² * τ^2 - 4.0im * τ * bath.ϵ)
    ohmic_part * slow_part
end

"""
    ohmic_correlation(τ, bath::HybridOhmicBath[, a=1])

Calculate the Ohmic part of polaron correlation function of HybridOhmicBath 'bath' at time 'τ' with relative strength `a`. The effective strength is `` a * η ``.
"""
function ohmic_correlation(τ, bath::HybridOhmicBath, a=1)
    η = a * bath.η
    res = (1+1.0im*bath.ωc*τ)^(-4*η)
    if !isapprox(τ, 0, atol = 1e-9)
        x = π * τ / bath.β
        res *= ( x / sinh(x) )^(4*η)
    end
    res
end

"""
    reorganization_energy(bath::HybridOhmicBath)

Calculate the total reorganization energy of HybridOhmicBath `bath`: ``ϵ = ϵ_L + ϵ_H``.
"""
function reorganization_energy(bath::HybridOhmicBath)
    bath.ϵ + 4 * bath.η * bath.ωc
end

function Base.show(io::IO, ::MIME"text/plain", m::HybridOhmicBath)
    print(io, "Hybrid Ohmic bath instance:\n", "W (mK): ", freq_2_temperature(m.W/2/pi), "\n",
        "ϵ (GHz): ", m.ϵ/2/pi ,"\n",
        "η (unitless): ", m.η, "\n", "ωc (GHz): ", m.ωc/pi/2, "\n", "T (mK): ", beta_2_temperature(m.β))
end

"""
    GH(ω, bath::HybridOhmicBath[, a=1])

High frequency noise spectrum of the HybridOhmicBath `bath` with relative strength `a`.
"""
function GH(ω, bath::HybridOhmicBath, a = 1)
    η = a * bath.η
    S0 = 8* pi* η / bath.β
    if isapprox(ω, 0, atol=1e-8)
        return 4 / S0
    else
        γ² = (S0 / 2)^2
        return 8 * pi * η * ω * exp(-abs(ω)/bath.ωc) / (1 - exp(-bath.β*ω)) / (ω^2 + γ²)
    end
end

"""
    GL(ω, bath::HybridOhmicBath[, a=1])

Low frequency noise specturm of the HybridOhmicBath `bath` with relative strength `a`.
"""
function GL(ω, bath::HybridOhmicBath, a=1)
    W² = a * bath.W^2
    ϵ = a * bath.ϵ
    sqrt(π/2/W²) * exp(-(ω-4*ϵ)^2/8/W²)
end

function convolution_rate(tf, sys, bath::HybridOhmicBath)
    Γ10 = []
    Γ01 = []
    ϵ = reorganization_energy(bath)
    for i in eachindex(sys.s)
        T_bar = sys.T[i] - 1.0im*sys.G[i] / tf - sys.d[i] * ϵ
        A = abs2(T_bar - sys.ω[i] * sys.c[i] / sys.a[i])
        B = (sys.a[i] * sys.b[i] - abs2(sys.c[i])) / sys.a[i]^2
        Δ²(ω) = A + B * (ω^2 + sys.a[i]*bath.W^2)
        integrand_12(ω) = Δ²(ω) * GL(sys.ω[i]-ω, bath, sys.a[i]) * GH(ω, bath, sys.a[i])
        integrand_21(ω) = Δ²(ω) * GL(-sys.ω[i]-ω, bath, sys.a[i]) * GH(ω, bath, sys.a[i])
        push!(Γ10, integrate_1d(integrand_12, -Inf, Inf)[1])
        push!(Γ01, integrate_1d(integrand_21, -Inf, Inf)[1])
    end
    Γ10/2/pi, Γ01/2/pi
end
