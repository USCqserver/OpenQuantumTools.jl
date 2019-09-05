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
    print(
        io,
        "Hybrid Ohmic bath instance:\n",
        "W (mK): ",
        freq_2_temperature(m.W / 2 / pi),
        "\n",
        "ϵl (GHz): ",
        m.ϵl / 2 / pi,
        "\n",
        "η (unitless): ",
        m.η / 2 / π,
        "\n",
        "ωc (GHz): ",
        m.ωc / pi / 2,
        "\n",
        "T (mK): ",
        beta_2_temperature(m.β),
        "\n",
        "ϵ (GHz): ",
        m.ϵ / pi / 2
    )
end

"""
    HybridOhmic(W, η, fc, T)

Construct HybridOhmicBath object with parameters in physical units. `W`: MRT width (mK); `η`: interaction strength (unitless); `fc`: Ohmic cutoff frequency (GHz); `T`: temperature (mK).
"""
function HybridOhmic(W, η, fc, T)
    # scaling of W to the unit of angular frequency
    W = 2 * pi * temperature_2_freq(W)
    # scaling of η because different definition of Ohmic spectrum
    η = 2 * π * η
    β = temperature_2_beta(T)
    ωc = 2 * π * fc
    ϵl = W^2 * β / 2
    ϵ = ϵl + η * ωc / 2 / π
    width_h = η / β / 2
    width_l = sqrt(2 * log(2)) * W
    HybridOhmicBath(W, ϵl, ϵ, η, ωc, β, width_h, width_l)
end


function MRT_Γ(ϵ, Δ, bath::HybridOhmicBath)
    integrand(τ) = p_correlation(τ, ϵ, bath)
    res, err = quadgk(integrand, -Inf, Inf, rtol=1e-8, atol=1e-8)
    Δ^2 * res, Δ^2 * err
end

function p_correlation(τ, ϵ, bath)
    η = bath.η / 2 / π
    ϵl = bath.ϵl
    W² = bath.W^2
    ohmic_part = (1 + 1.0im * bath.ωc * τ)^(-η)
    if !isapprox(τ, 0, atol = 1e-9)
        x = π * τ / bath.β
        ohmic_part *= (x / sinh(x))^(η)
    end
    slow_part = exp(-W² * τ^2 / 2 - 1.0im * τ * (ϵl-ϵ))
    ohmic_part * slow_part
end


"""
    polaron_correlation(τ, bath::HybridOhmicBath, a=1)

Calculate polaron transformed correlation function of HybridOhmicBath 'bath' at time 'τ' with relative strength `a`. The effective strength will be `` a * W^2`` and `` a * η ``.
"""
function polaron_correlation(τ, bath::HybridOhmicBath, a = 1)
    η = a * bath.η / 2 / π
    ϵ = a * bath.ϵl
    W² = a * bath.W^2
    ohmic_part = (1 + 1.0im * bath.ωc * τ)^(-η)
    if !isapprox(τ, 0, atol = 1e-9)
        x = π * τ / bath.β
        ohmic_part *= (x / sinh(x))^(η)
    end
    slow_part = exp(-W² * τ^2 / 2 - 1.0im * τ * ϵ)
    ohmic_part * slow_part
end

"""
    Gₕ(ω, bath::HybridOhmicBath, a=1)

High frequency noise spectrum of the HybridOhmicBath `bath` with relative strength `a`.
"""
function Gₕ(ω, bath::HybridOhmicBath, a = 1)
    η = a * bath.η
    S0 = η / bath.β
    if isapprox(ω, 0, atol = 1e-8)
        return 4 / S0
    else
        γ² = (S0 / 2)^2
        return η * ω * exp(-abs(ω) / bath.ωc) / (1 - exp(-bath.β * ω)) /
               (ω^2 + γ²)
    end
end

"""
    Gₗ(ω, bath::HybridOhmicBath, a=1)

Low frequency noise specturm of the HybridOhmicBath `bath` with relative strength `a`.
"""
function Gₗ(ω, bath::HybridOhmicBath, a = 1)
    W² = a * bath.W^2
    ϵ = a * bath.ϵl
    sqrt(2 * π / W²) * exp(-(ω - ϵ)^2 / 2 / W²)
end

"""
    bloch_rate(i, tf, sys, bath::HybridOhmicBath)

Calculate the relaxation rate of polaron transformed ME in the Bloch-Redfield limit.
"""
function bloch_rate(i, tf, sys, bath::HybridOhmicBath)
    ω₀₁ = sys.ω[i] - sys.a[i] * bath.ϵl
    ω₁₀ = -sys.ω[i] - sys.a[i] * bath.ϵl
    W² = bath.W^2
    ω²₀₁ = ω₀₁^2
    ω²₁₀ = ω₁₀^2
    a = sys.a[i]
    b = sys.b[i]
    c = sys.c[i]
    T̃ = sys.T[i] - 1.0im * sys.G[i] / tf - sys.d[i] * bath.ϵ
    S₀₁ = Sₕ(ω₀₁, bath)
    S₁₀ = Sₕ(ω₁₀, bath)
    temp = (a * (abs2(T̃) + b * W²) + b * ω²₀₁ - 2 * ω₀₁ * real(T̃) * c -
            abs2(c) * W²)
    Γ₁₀ = temp * S₀₁ / (ω²₀₁ + a^2 * bath.η^2 / bath.β^2 / 4)
    temp = (a * (abs2(T̃) + b * W²) + b * ω²₁₀ - 2 * ω₁₀ * real(T̃) * c -
            abs2(c) * W²)
    Γ₀₁ = temp * S₁₀ / (ω²₁₀ + a^2 * bath.η^2 / bath.β^2 / 4)
    Γ₀₁, Γ₁₀
end

"""
    direct_integrate(i, tf, sys, bath::HybridOhmicBath)

Calculate the relaxation rate of polaron transformed ME by directly integrating the convolution formula.
"""
function direct_integrate(i, tf, sys, bath::HybridOhmicBath)
    # get the spectrum widths and centers
    # we use 3σ for Gaussian profile
    μ₀₁ = sys.ω[i] - sys.a[i] * bath.ϵl
    μ₁₀ = -sys.ω[i] - sys.a[i] * bath.ϵl
    σ = sqrt(sys.a[i]) * bath.width_l
    γ = sys.a[i] * bath.width_h
    integration_range_01 = sort([
        μ₀₁ - 3 * σ,
        μ₀₁,
        μ₀₁ + 3 * σ,
        -3 * γ,
        0,
        3 * γ,
        2 * σ,
        -2 * σ
    ])
    integration_range_10 = sort([
        μ₁₀ - 3 * σ,
        μ₁₀,
        μ₁₀ + 3 * σ,
        -3 * γ,
        0,
        3 * γ,
        2 * σ,
        -2 * σ
    ])
    #
    T̃ = sys.T[i] - 1.0im * sys.G[i] / tf - sys.d[i] * bath.ϵ
    A₀₁ = abs2(T̃ - sys.ω[i] * sys.c[i] / sys.a[i])
    A₁₀ = abs2(T̃ + sys.ω[i] * sys.c[i] / sys.a[i])
    B = (sys.a[i] * sys.b[i] - abs2(sys.c[i])) / sys.a[i]^2
    Δ²₀₁(ω) = A₀₁ + B * (ω^2 + sys.a[i] * bath.W^2)
    Δ²₁₀(ω) = A₁₀ + B * (ω^2 + sys.a[i] * bath.W^2)
    integrand_01(ω) =
        Δ²₀₁(ω) * Gₗ(sys.ω[i] - ω, bath, sys.a[i]) * Gₕ(ω, bath, sys.a[i])
    integrand_10(ω) =
        Δ²₁₀(ω) * Gₗ(-sys.ω[i] - ω, bath, sys.a[i]) * Gₕ(ω, bath, sys.a[i])
    Γ₁₀, err₁₀ = quadgk(
        integrand_01,
        integration_range_01...,
        rtol = 1e-6,
        atol = 1e-6
    )
    Γ₀₁, err₀₁ = quadgk(
        integrand_10,
        integration_range_10...,
        rtol = 1e-6,
        atol = 1e-6
    )
    (Γ₀₁ / 2 / π, Γ₁₀ / 2 / π), (err₁₀, err₀₁)
end

function integration_range(i, integrand_01, integrand_10, sys, bath)
    μ₀₁ = sys.ω[i] - sys.a[i] * bath.ϵl
    μ₁₀ = -sys.ω[i] - sys.a[i] * bath.ϵl
    σ = sqrt(sys.a[i]) * bath.width_l
    γ = sys.a[i] * bath.width_h
    integration_range_01 = sort([
        μ₀₁ - 3 * σ,
        μ₀₁,
        μ₀₁ + 3 * σ,
        -3 * γ,
        0,
        3 * γ,
        4 * γ
    ])
    integration_range_10 = sort([
        μ₁₀ - 3 * σ,
        μ₁₀,
        μ₁₀ + 3 * σ,
        -3 * γ,
        0,
        3 * γ,
        4 * γ
    ])
    # range for 01 term
    integrand_01(μ₀₁)
    # range for 10 term
    integrand_10()
end

function spectrum_info(bath::HybridOhmicBath, a = 1)
    Dict(
        "low_freq_center" => a * bath.ϵl,
        "low_freq_height" => sqrt(2 * π / a) / bath.W,
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

"""
    Sₕ(ω, bath::HybridOhmicBath)

The corresponding Ohmic spectrum of `HybridOhmicBath`.
"""
function Sₕ(ω, bath::HybridOhmicBath)
    isapprox(ω, 0, atol = 1e-8) ? bath.η / bath.β :
    bath.η * ω * exp(-abs(ω) / bath.ωc) / (1 - exp(-bath.β * ω))
end
