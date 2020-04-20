struct Interaction
    coupling::AbstractCouplings
    bath::AbstractBath
end

function create_redfield(inter::Interaction, unitary, tf)
    create_redfield(inter.coupling, unitary, tf, inter.bath)
end

function build_davies(inter::Interaction, unitary, tf)
    build_davies(inter.coupling, inter.bath, ω_hint, lambshift)
end

struct InteractionSet{T<:Tuple}
    interactions::T
end

function InteractionSet(inters::Interaction...)
    InteractionSet(inters)
end

function create_redfield(inter::InteractionSet, unitary, tf)
    RedfieldSet([
        create_redfield(i.coupling, unitary, tf, i.bath)
        for i in inter.interactions
    ]...)
end


function build_davies(inter::InteractionSet, ω_range, lambshift)
    error("InteractionSet is not supported for adiabatic master equation solver.")
end


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
