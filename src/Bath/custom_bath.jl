mutable struct CustomBath <: AbstractBath
    cfun
    γ
end


function CustomBath(; correlation = nothing, spectrum = nothing)
    CustomBath(correlation, spectrum)
end


function correlation(τ, bath::CustomBath)
    bath.cfun(τ)
end


function create_redfield(
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


function create_redfield(
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
