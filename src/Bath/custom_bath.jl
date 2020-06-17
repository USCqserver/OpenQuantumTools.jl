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

function build_redfield(
    coupling,
    unitary,
    tf::Real,
    bath::CustomBath;
    atol = 1e-8,
    rtol = 1e-6,
)
    if bath.cfun == nothing
        error("Correlation function is not defined for the bath.")
    end
    cfun(s) = bath.cfun(s * tf)
    Redfield(coupling, unitary, cfun, atol = atol, rtol = rtol)
end

function build_redfield(
    coupling,
    unitary,
    tf::UnitTime,
    bath::CustomBath;
    atol = 1e-8,
    rtol = 1e-6,
)
    if bath.cfun == nothing
        error("Correlation function is not defined for the bath.")
    end
    cfun(t) = bath.cfun(t)
    Redfield(coupling, unitary, cfun, atol = atol, rtol = rtol)
end

function build_davies(coupling, bath::CustomBath, ω_range, lambshift)
    if bath.γ == nothing
        error("Noise spectrum is not defined for the bath.")
    end
    if lambshift == true
        if isempty(ω_range)
            S_loc = (ω) -> lambshift(ω, bath.γ)
        else
            s_list = [S(ω, bath) for ω in ω_range]
            S_loc = construct_interpolations(ω_range, s_list)
        end
    else
        S_loc = (ω) -> 0.0
    end
    DaviesGenerator(coupling, bath.γ, S_loc)
end
