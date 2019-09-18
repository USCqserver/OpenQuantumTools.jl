"""
$(TYPEDEF)

A hybrid noise model with both low and high frequency noise. The high frequency noise is characterized by Ohmic bath and the low frequence noise is characterized by the MRT width `W`.

$(FIELDS)
"""
struct HybridOhmicBath
    """MRT width (2π GHz)"""
    W::Float64
    """low spectrum reorganization energy (2π GHz)"""
    ϵl::Float64
    """total reorganization energy (2π GHz)"""
    ϵ::Float64
    """strength of high frequency Ohmic bath"""
    η::Float64
    """cutoff frequency"""
    ωc::Float64
    """inverse temperature"""
    β::Float64
    """half width at half maximum for high frequency Ohmic bath"""
    width_h::Float64
    """half width at half maximum for low frequency bath"""
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
    Gl(ω) = Gₗ(ω, bath)
    Gh(ω) = Gₕ(ω, bath)
    integrand(ω) = Gl(ϵ-ω)*Gh(ω)
    a, b = sort([0.0, ϵ])
    res, err = quadgk(integrand, -Inf, a, b, Inf)
    Δ^2*res/8/π, Δ^2*err/8/π
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
