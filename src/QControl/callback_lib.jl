"""This is a modified version of PresetTimeCallback object in DiffEqCallbacks.jl"""
function PresetTimeCallback(tstops,user_affect!;
                            initialize = INITIALIZE_DEFAULT, kwargs...)
    condition = function (u, t, integrator)
      t in scale_tstops(integrator.p.tf, tstops)
    end

    # Call f, update tnext, and make sure we stop at the new tnext
    affect! = function (integrator)
        user_affect!(integrator)
        nothing
    end

    # Initialization: first call to `f` should be *before* any time steps have been taken:
    initialize_preset = function (c, u, t, integrator)
        initialize(c, u, t, integrator)
        scaled_tstops = scale_tstops(integrator.p.tf, tstops)
        add_tstop!.((integrator,), scaled_tstops)
        if t ∈ scaled_tstops
            user_affect!(integrator)
        end
    end
    DiscreteCallback(condition, affect!; initialize = initialize_preset, kwargs...)
end

scale_tstops(tf::UnitTime, tstops) = tf * tstops
scale_tstops(tf::Real, tstops) = tstops


function empty_affect!(integrator)
    u_modified!(integrator,false)
end
