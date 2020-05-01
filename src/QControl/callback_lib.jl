"""This is a modified version of PresetTimeCallback object in DiffEqCallbacks.jl"""
function PresetTimeCallback(
    tstops,
    user_affect!;
    initialize = INITIALIZE_DEFAULT,
    filter_tstops = true,
    kwargs...,
)
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
        if filter_tstops
            tdir = integrator.tdir
            _tstops = scaled_tstops[@.(
                (tdir * tstops >= tdir * integrator.sol.prob.tspan[1]) *
                (tdir * tstops < tdir * integrator.sol.prob.tspan[2])
            )]
            add_tstop!.((integrator,), _tstops)
        else
            add_tstop!.((integrator,), scaled_tstops)
        end
        if t ∈ scaled_tstops
            user_affect!(integrator)
        end
    end
    DiscreteCallback(
        condition,
        affect!;
        initialize = initialize_preset,
        kwargs...,
    )
end

scale_tstops(tf::UnitTime, tstops) = tf * tstops
scale_tstops(tf::Real, tstops) = tstops


function positivity_check_affect(u, t, integrator)
    if !check_positivity(u)
        @warn "The density matrix becomes negative at time $t."
        terminate!(integrator)
    end
end


function positivity_check_affect(
    u::AbstractArray{T,1},
    t,
    integrator,
) where {T<:Number}
    dim = convert(Int, sqrt(length(u)))
    if !check_positivity(reshape(u, dim, dim))
        @warn "The density matrix becomes negative at time $t."
        terminate!(integrator)
    end
end


"""This is a modified version of IterativeCallback object in DiffEqCallbacks.jl"""
function IterativeCallback(
    time_choice,
    user_affect!,
    tType = Float64;
    initialize = INITIALIZE_DEFAULT,
    initial_affect = false,
    kwargs...,
)
    # Value of `t` at which `f` should be called next:
    tnext = Ref{Union{Nothing,eltype(tType)}}(typemax(tType))
    condition = function (u, t, integrator)
        t == tnext[]
    end

    # Call f, update tnext, and make sure we stop at the new tnext
    affect! = function (integrator)
        user_affect!(integrator)

        # Schedule next call to `f` using `add_tstops!`, but be careful not to keep integrating forever
        tnew = time_choice(integrator)
        tnew === nothing && (tnext[] = tnew; return)
        tstops = integrator.opts.tstops
        for i = length(tstops):-1:1 # reverse iterate to encounter large elements earlier
            #=
            Okay yeah, this is nasty
            the comparer is always less than for type stability, so in order
            for this to actually check the correct direction we multiply by
            tdir
            =#
            if compare(
                tstops.comparer,
                integrator.tdir * tnew,
                integrator.tdir * tstops.valtree[i],
            ) # TODO: relying on implementation details
                tnext[] = tnew
                add_tstop!(integrator, tnew)
                break
            elseif tstops.valtree[i] == tnew
                # If it's already a tstop, no need to re-add! This is for the final point
                tnext[] = tnew
            end
        end
        nothing
    end

    # Initialization: first call to `f` should be *before* any time steps have been taken:
    initialize_iterative = function (c, u, t, integrator)
        initialize(c, u, t, integrator)
        if initial_affect
            tnext[] = t
            affect!(integrator)
        else
            tnext[] = time_choice(integrator)
            if tnext[] != nothing
                add_tstop!(integrator, tnext[])
            end
        end
    end
    DiscreteCallback(
        condition,
        affect!;
        initialize = initialize_iterative,
        kwargs...,
    )
end


function AMEJumpCallback()
    condition = function (u, t, integrator)
        real(u.x' * u.x) - u.r
    end

    affect! = function (integrator)
        A = QTBase.ame_jump(
            integrator.p.control,
            integrator.u,
            integrator.p.tf,
            integrator.t,
        )
        new_r = rand()
        new_state = normalize(A * integrator.u)
        for c in full_cache(integrator)
            c.x .= new_state
            c.r = new_r
        end
    end

    ContinuousCallback(condition, affect!)
end
