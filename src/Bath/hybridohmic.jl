"""
    HybridOhmicBath

A hybrid noise model with both low and high frequency noise. The high frequency noise is characterized by Ohmic bath and the low frequence noise is characterized by the MRT width `W`.

**Fields**
- `W` -- MRT width (2π GHz)
- `ϵl` -- low spectrum reorganization energy (2π GHz)
- `ϵ` -- total reorganization energy (2π GHz)
- `η` -- strength of high frequency Ohmic bath
- `ωc` -- cutoff frequence
- `β` -- inverse temperature
- `width_h` -- half width at half maximum for high frequency Ohmic bath
- `width_l` -- half width at half maximu for low frequency slow bath
"""
struct HybridOhmicBath
    W::Float64
    ϵl::Float64
    ϵ::Float64
    η::Float64
    ωc::Float64
    β::Float64
    width_h::Float64
    width_l::Float64
end

function Base.show(io::IO, ::MIME"text/plain", m::HybridOhmicBath)
    print(io, "Hybrid Ohmic bath instance:\n", "W (mK): ", freq_2_temperature(m.W/4/pi), "\n",
        "ϵl (GHz): ", m.ϵl/8/pi ,"\n",
        "η (unitless): ", m.η/8/π, "\n", "ωc (GHz): ", m.ωc/pi/2, "\n", "T (mK): ", beta_2_temperature(m.β), "\n", "ϵ (GHz): ", m.ϵ/pi/8)
end

"""
    HybridOhmic(W, η, fc, T)

Construct HybridOhmicBath object with parameters in physical units. `W`: MRT width (mK); `η`: interaction strength (unitless); `fc`: Ohmic cutoff frequency (GHz); `T`: temperature (mK).
"""
function HybridOhmic(W, η, fc, T)
    # scaling of W comes from the definition of interaction term
    W = 2 * 2 * pi * temperature_2_freq(W)
    # scaling of η comes from the definition of interaction term and the ohmic spectrum
    η = 4 * 2 * π * η

    β = temperature_2_beta(T)
    ωc = 2 * π * fc
    ϵl = W^2 * β / 2
    ϵ = ϵl + η * ωc / 2 / π
    width_h = η / β / 2
    width_l = sqrt(2 * log(2)) * W
    HybridOhmicBath(W, ϵl, ϵ, η, ωc, β, width_h, width_l)
end

"""
    polaron_correlation(τ, bath::HybridOhmicBath, a=1)

Calculate polaron transformed correlation function of HybridOhmicBath 'bath' at time 'τ' with relative strength `a`. The effective strength will be `` a * W^2`` and `` a * η ``.
"""
function polaron_correlation(τ, bath::HybridOhmicBath, a=1)
    η = a * bath.η
    ϵ = a * bath.ϵl
    W² = a * bath.W^2
    ohmic_part = (1+1.0im*bath.ωc*τ)^(-η)
    if !isapprox(τ, 0, atol = 1e-9)
        x = π * τ / bath.β
        ohmic_part *= ( x / sinh(x) )^(η)
    end
    slow_part = exp( - W² * τ^2 - 1.0im * τ * ϵ)
    ohmic_part * slow_part
end

"""
    GH(ω, bath::HybridOhmicBath, a=1)

High frequency noise spectrum of the HybridOhmicBath `bath` with relative strength `a`.
"""
function GH(ω, bath::HybridOhmicBath, a = 1)
    η = a * bath.η
    S0 = η / bath.β
    if isapprox(ω, 0, atol=1e-8)
        return 4 / S0
    else
        γ² = (S0 / 2)^2
        return η * ω * exp(-abs(ω)/bath.ωc) / (1 - exp(-bath.β*ω)) / (ω^2 + γ²)
    end
end

"""
    GL(ω, bath::HybridOhmicBath, a=1)

Low frequency noise specturm of the HybridOhmicBath `bath` with relative strength `a`.
"""
function GL(ω, bath::HybridOhmicBath, a=1)
    W² = a * bath.W^2
    ϵ = a * bath.ϵl
    sqrt(2*π/W²) * exp(-(ω-ϵ)^2/2/W²)
end

function spectrum_info(bath::HybridOhmicBath, a = 1)
    Dict(
        "low_freq_center" => a * bath.ϵl,
        "low_freq_height" => sqrt(2*π/a)/bath.W,
        "low_freq_FWHM" => bath.width_l * sqrt(a) * 2,
        "high_freq_height" => 4 * bath.β / a / bath.η,
        "high_freq_FWHM" => bath.width_h * a * 2
    )
end

"""
    half_width_half_maximum(a, bath::HybridOhmicBath)

Calculate the half width at half maximum of both the low/high frequency spectrum function ``G_L(ω)`` and ``G_H(ω)``.
"""
function half_width_half_maximum(a, bath::HybridOhmicBath)
    Wh = bath.width_h * a
    Wl = bath.width_l * sqrt(a)
    Wh, Wl
end

function tunneling_Δ(ω, i, sys, tf, bath::HybridOhmicBath)
    T_bar = sys.T[i] - 1.0im*sys.G[i] / tf - sys.d[i] * bath.ϵ
    A = abs2(T_bar - sys.ω[i] * sys.c[i] / sys.a[i])
    B = (sys.a[i] * sys.b[i] - abs2(sys.c[i])) / sys.a[i]^2
    A + B * (ω^2 + sys.a[i]*bath.W^2)
end

"""
    bloch_rate(i, tf, sys, bath::HybridOhmicBath)

Calculate the relaxation rate of polaron transformed ME in the Bloch-Redfield limit.
"""
function bloch_rate(i, tf, sys, bath::HybridOhmicBath)
    ω = sys.ω[i] - sys.a[i] * bath.ϵl
    W² = bath.W^2
    ω² = ω^2
    a = sys.a[i]
    b = sys.b[i]
    c = sys.c[i]
    T̃ = sys.T[i] - 1.0im * sys.G[i] / tf - sys.d[i] * bath.ϵ
    SH = bath.η * ω * exp(-abs(ω)/bath.ωc) / (1 - exp(-bath.β*ω))
    temp1 = (a * (T̃^2 + b*W²) + b*ω² - 2*ω*T̃*c - abs2(c)*W²)
    temp1 * SH / (ω² + a^2*bath.η^2/bath.β^2/4)
end

function direct_integrate(i, tf, sys, bath::HybridOhmicBath)
    T̃ = sys.T[i] - 1.0im * sys.G[i] / tf - sys.d[i] * bath.ϵ
    A = abs2(T̃ - sys.ω[i] * sys.c[i] / sys.a[i])
    B = (sys.a[i] * sys.b[i] - abs2(sys.c[i])) / sys.a[i]^2
    Δ²(ω) = A + B * (ω^2 + sys.a[i]*bath.W^2)
    low_width = bath.width_l * sqrt(a) * 2
    high_width = bath.width_h * a * 2
end

function convolution_rate(tf, sys, bath::HybridOhmicBath)
    Γ10 = []
    Γ01 = []
    ϵ = bath.ϵ
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
