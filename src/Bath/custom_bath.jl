"""
$(TYPEDEF)

An custum bath object defined by the two-point correlation function and the corresponding spectrum.

$(FIELDS)
"""
mutable struct CustomBath <: AbstractBath
    """correlation function"""
    cfun::Any
    """spectrum"""
    γ::Any
end

CustomBath(; correlation = nothing, spectrum = nothing) =
    CustomBath(correlation, spectrum)
correlation(τ, bath::CustomBath) = bath.cfun(τ)
spectrum(ω, bath::CustomBath) = bath.γ(ω)
γ(ω, bath::CustomBath) = bath.γ(ω)
S(w, bath::CustomBath; atol = 1e-7) =
    lambshift(w, (ω) -> spectrum(ω, bath), atol = atol)
build_correlation(bath::CustomBath, tf::Real) =
    bath.cfun == nothing ? error("Correlation function is not specified.") :
    (s) -> bath.cfun(s * tf)
build_correlation(bath::CustomBath, ::UnitTime) =
    bath.cfun == nothing ? error("Correlation function is not specified.") :
    bath.cfun
build_spectrum(bath::CustomBath) =
    bath.γ == nothing ? error("Noise spectrum is not specified.") : bath.γ
