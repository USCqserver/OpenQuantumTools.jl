"""
$(TYPEDEF)

An object to hold system operator and the corresponding bath object.

$(FIELDS)
"""
struct Interaction
    """system operator"""
    coupling::AbstractCouplings
    """bath coupling to the system operator"""
    bath::AbstractBath
end


function build_redfield(
    inter::Interaction,
    unitary,
    tf;
    atol = 1e-8,
    rtol = 1e-6,
)
    build_redfield(
        inter.coupling,
        unitary,
        tf,
        inter.bath,
        atol = atol,
        rtol = rtol,
    )
end


build_davies(inter::Interaction, ω_hint, lambshift) =
    build_davies(inter.coupling, inter.bath, ω_hint, lambshift)


"""
$(TYPEDEF)

An container for different system-bath interactions.

$(FIELDS)
"""
struct InteractionSet{T<:Tuple}
    """A tuple of Interaction"""
    interactions::T
end


function InteractionSet(inters::Interaction...)
    InteractionSet(inters)
end


function build_redfield(
    inter::InteractionSet,
    unitary,
    tf;
    atol = 1e-8,
    rtol = 1e-6,
)
    RedfieldSet([
        build_redfield(i, unitary, tf, atol = atol, rtol = rtol)
        for i in inter.interactions
    ]...)
end


build_davies(inter::InteractionSet, ω_range, lambshift) =
    error("InteractionSet is not supported for adiabatic master equation solver.")


function build_ame_trajectory_control_from_interactions(
    inter::InteractionSet,
    ω_hint,
    lambshift,
    lvl,
    tf,
    H,
    ame_trajectory_de_field,
    fluctuator_de_field,
)
    a_control = []
    f_control = []
    opensys = nothing
    for i in inter.interactions
        if typeof(i.bath) <: EnsembleFluctuator
            c = FluctuatorControl(
                tf,
                length(i.coupling),
                i.bath,
                fluctuator_de_field,
            )
            push!(f_control, c)
            opensys = StochasticNoise(i.coupling, fluctuator_de_field)
        else
            d = build_davies(i.coupling, i.bath, ω_hint, lambshift)
            op = AMETrajectoryOperator(H, d, lvl)
            push!(a_control, AMETrajectoryControl(op, ame_trajectory_de_field))
        end
    end
    if length(f_control) > 1
        error("Only single fluctuator ensemble is supported.")
    end
    if length(a_control) > 1
        error("Multi-axis bath is not yet supported.")
    end
    ControlSet(a_control..., f_control...), opensys
end
