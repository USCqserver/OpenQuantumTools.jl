function build_ss_control_from_interactions(
    i::Interaction,
    tf,
    field::Union{Nothing,Symbol},
)
    control = FluctuatorControl(tf, length(i.coupling), i.bath, field)
    opensys = StochasticNoise(i.coupling, field)
    control, opensys
end

function build_ss_control_from_interactions(
    iset::InteractionSet,
    tf,
    field::Union{Nothing,Symbol},
)
    if length(iset) == 1
        build_ss_control_from_interactions(iset[1], tf, field)
    else
        throw(ArgumentError("Multiple interactions is not yet supported for stochastic Schrodinger equation solver."))
    end
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
            d = build_davies(i, ω_hint, lambshift)
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

function build_hybrid_redfield_control_from_interactions(
    inter::InteractionSet,
    unitary,
    tf,
    Ta,
    atol,
    rtol,
    fluctuator_de_field,
)
    f_control = []
    f_opensys = []
    r_opensys = []
    for i in inter.interactions
        if typeof(i.bath) <: EnsembleFluctuator
            c = FluctuatorControl(
                tf,
                length(i.coupling),
                i.bath,
                fluctuator_de_field,
            )
            push!(f_control, c)
            push!(f_opensys, StochasticNoise(i.coupling, fluctuator_de_field))
        else
            s = QTBase.redfield_from_interactions(
                i,
                unitary,
                tf,
                Ta,
                atol = atol,
                rtol = rtol,
            )
            push!(r_opensys, s)
        end
    end
    if isempty(f_control)
        error("No fluctuator ensemble detected. Use solve_redfield instead.")
    elseif length(f_control) > 1
        error("Only single fluctuator ensemble is supported.")
    end

    if length(r_opensys) == 1
        r_opensys = r_opensys[1]
    else
        r_opensys = RedfieldSet(r_opensys...)
    end
    # currently only sinlge fluctuator ensemble is supported
    f_control[1], f_opensys[1], r_opensys
end
