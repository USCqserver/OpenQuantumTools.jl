mutable struct CustomBath <: AbstractBath
    cfun
    γ
end


function CustomBath(;correlation=nothing, spectrum=nothing)
    CustomBath(correlation, spectrum)
end


function create_redfield(coupling, unitary, tf::Real, bath::CustomBath)
    if bath.cfun == nothing
        error("Correlation function is not defined for the bath.")
    end
    cfun(s) = bath.cfun(s * tf)
    Redfield(coupling, unitary, cfun)
end


function create_redfield(coupling, unitary, tf::UnitTime, bath::CustomBath)
    if bath.cfun == nothing
        error("Correlation function is not defined for the bath.")
    end
    cp(t) = coupling(t/tf)
    cfun(t) = bath.cfun(t)
    Redfield(cp, unitary, cfun)
end
