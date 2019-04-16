function Ohmic(η, fc, T)
    ωc = 2 * π * fc
    β = temperature_2_beta(T)
    OhmicBath(η, ωc, β)
end

function HybridOhmic(W, η, fc, T)
    W = 2 * pi * temperature_2_freq(W)
    β = temperature_2_beta(T)
    ϵ = W^2 * β / 2
    ωc = 2 * π * fc
    HybridOhmicBath(W, ϵ, η, ωc, β)
end

function Base.show(io::IO, ::MIME"text/plain", m::OhmicBath)
    print(io, "Ohmic bath instance:\n", "η (unitless): ", m.η, "\n", "ωc (GHz): ", m.ωc/pi/2,
        "\n", "T (mK): ", beta_2_temperature(m.β))
end

function Base.show(io::IO, ::MIME"text/plain", m::HybridOhmicBath)
    print(io, "Hybrid Ohmic bath instance:\n", "W (mK): ", freq_2_temperature(m.W/2/pi), "\n",
        "ϵ (GHz): ", m.ϵ/2/pi ,"\n",
        "η (unitless): ", m.η, "\n", "ωc (GHz): ", m.ωc/pi/2, "\n", "T (mK): ", beta_2_temperature(m.β))
end
