function DEFAULT_INITIALIZER end

"""
$(SIGNATURES)

Builds the callback to check positivity of the system density matrix during the evolution.
"""
function PositivityCheckCallback()
    affect! = function (u, t, integrator)
        if ndims(u) == 1
            dim = convert(Int, sqrt(length(u)))
            ρ = reshape(u, dim, dim)
        else
            ρ = u
        end
        if !check_positivity(ρ)
            @warn "The density matrix becomes negative at time $t."
            terminate!(integrator)
        end
        u_modified!(integrator, false)
    end
    FunctionCallingCallback(affect!, func_everystep=true, func_start=false)
end

"""
$(SIGNATURES)

Builds the callback to apply instantaneous pulses during the evolution.

...
# Arguments
- `tstops`: time points to apply pulses.
- `pulse_update`: udpate function for the state of the system. It taks two argument `pulse_update(c, i)` where `c` is the state vector or density matrix and `i` is the index of the pulse being applied.
...
"""
function InstPulseCallback(tstops, pulse_update)
    state = Ref{Int}(1)
    affect! = function (integrator)
        for c in full_cache(integrator)
            pulse_update(c, state[])
        end
        state[] += 1
    end

    initialize = function (c, u, t, integrator)
        state[] = 1
        u_modified!(integrator, false)
    end
    PresetTimeCallback(tstops, affect!, initialize=initialize)
end

function FluctuatorCallback(F::OpenQuantumBase.FluctuatorLiouvillian, initialize, save_positions)
    time_choice = function (integrator)
        next_t = integrator.t + F.next_τ
        if next_t > integrator.sol.prob.tspan[2]
            return nothing
        else
            return next_t
        end
    end

    affect! = function (integrator)
        F.n = sum(F.b0, dims=1)[:]
        OpenQuantumBase.next_state!(F)
        u_modified!(integrator, false)
    end

    initialize_fluc = function (cb, u, t, integrator)
        OpenQuantumBase.reset!(F, initialize)
        u_modified!(integrator, false)
    end

    IterativeCallback(
        time_choice,
        affect!,
        initialize=initialize_fluc,
        save_positions=save_positions,
    )
end

function LindbladJumpCallback(save_positions)
    r = Ref{Float64}(rand(Float64))
    condition = function (u, t, integrator)
        real(u' * u) - r[]
    end

    affect! = function (integrator)
        A = OpenQuantumBase.lindblad_jump(
            integrator.p.L,
            integrator.u,
            integrator.p,
            integrator.t,
        )
        r[] = rand()
        new_state = normalize(A * integrator.u)
        for c in full_cache(integrator)
            c .= new_state
        end
    end

    initialize = function (c, u, t, integrator)
        r[] = rand()
        u_modified!(integrator, false)
    end

    ContinuousCallback(condition, affect!, save_positions=save_positions, initialize=initialize)
end

"""
$(SIGNATURES)

Create callback to apply manifold retraction after each step of the ODE solver.
"""
function ManifoldRetractionCallback()
    affect! = function (integrator)
        rmul!(integrator.u, 1 / norm(integrator.u))
    end
    condition = (u, t, integrator) -> true
    DiscreteCallback(condition, affect!, save_positions=(false, save))
end
