"""
    function solve_ame(
        A::Annealing,
        tf::Real;
        dimensionless_time::Bool = true,
        ω_hint = [],
        lambshift::Bool = true,
        lvl::Int = size(A.H, 1),
        tstops = Float64[],
        vectorize::Bool = false,
        de_array_constructor = nothing,
        kwargs...,
    )

Solve the adiabatic master equation for `Annealing` defined by `A` with total annealing time `tf`.

...
# Arguments
- `A::Annealing`: the Annealing object.
- `tf::Real`: the total annealing time.
- `dimensionless_time::Bool=true`: flag variable which, when set to true, informs the solver to work with dimensionless time.
- `ω_hint=[]` : grid for precalculating the lambshift; skip the precalculation if empty.
- 'lambshift::Bool=true' : whether to include Lambshift in the calculation.
- `lvl::Int=size(A.H, 1)` : number of levels to keep. The default value is the dimension for the Hamiltonian.
- `tstops` : extra times that the timestepping algorithm must step to.
- `vectorize::Bool = false`: whether to vectorize the density matrix.
- `de_array_constructor = nothing`: the converting function if using `DEDataArray` type.
- `kwargs` : other keyword arguments supported by DifferentialEquations.jl
...
"""
function solve_ame(
    A::Annealing,
    tf::Real;
    dimensionless_time::Bool = true,
    ω_hint = [],
    lambshift::Bool = true,
    lvl::Int = size(A.H, 1),
    tstops = Float64[],
    vectorize::Bool = false,
    de_array_constructor = nothing,
    kwargs...,
)
    tf, u0, tstops = __init(
        A,
        tf,
        dimensionless_time,
        :m,
        tstops,
        de_array_constructor,
        vectorize = vectorize,
    )
    davies = build_davies(A.interactions, ω_hint, lambshift)
    if vectorize
        error("Vectorization is not yet supported for adiabatic master equation.")
    end
    f = AMEDiffEqOperator(A.H, davies, lvl)
    reset!(A.control)
    callback = ame_build_callback(A.control)
    p = ODEParams(tf; control = A.control)
    prob = ODEProblem(f, u0, (p) -> scaling_tspan(p.tf, A.sspan), p)
    solve(
        prob;
        alg_hints = [:nonstiff],
        callback = callback,
        tstops = tstops,
        kwargs...,
    )
end


ame_build_callback(control) = nothing
ame_build_callback(control::Union{InstPulseControl,InstDEPulseControl}) =
    build_callback(control, pulse_on_density!)

function build_ensemble_problem_ame_trajectory(
    A::Annealing,
    tf::Real,
    output_func,
    reduction;
    dimensionless_time::Bool = true,
    ω_hint = [],
    lambshift::Bool = true,
    lvl::Int = size(A.H, 1),
    de_array_constructor = nothing,
    ame_trajectory_de_field = nothing,
    fluctuator_de_field = nothing,
    initializer = DEFAULT_INITIALIZER,
    tstops = Float64[],
    kwargs...,
)
    # Put field names for ame_trajectory `r` and fluctuator `n` into `additional_symbol`
    # TODO: Do we really need to use keyword argument for those field names
    additional_symbol = []
    if ame_trajectory_de_field != nothing
        push!(additional_symbol, ame_trajectory_de_field)
    end
    if fluctuator_de_field != nothing
        push!(additional_symbol, fluctuator_de_field)
    end

    tf, u0, tstops = __init(
        A,
        tf,
        dimensionless_time,
        :v,
        tstops,
        de_array_constructor,
        needed_symbol = additional_symbol,
    )

    control, opensys = build_ame_trajectory_control_from_interactions(
        A.interactions,
        ω_hint,
        lambshift,
        lvl,
        tf,
        A.H,
        ame_trajectory_de_field,
        fluctuator_de_field,
    )
    control = A.control == nothing ? control : ControlSet(control, A.control)
    p = ODEParams(tf; control = control, opensys = opensys)
    callback = ame_trajectory_build_callback(control)
    ff = ame_trajectory_build_ode_function(A.H, control)
    prob = ODEProblem{true}(ff, u0, (p) -> scaling_tspan(p.tf, A.sspan), p)

    prob_func = build_prob_func(initializer)

    ensemble_prob = EnsembleProblem(
        prob;
        prob_func = prob_func,
        output_func = output_func,
        reduction = reduction,
    )

    ensemble_prob, callback, tstops
end


function ame_trajectory_build_callback(control::ControlSet)
    callbacks = [
        ame_trajectory_build_callback(v, k)
        for (k, v) in zip(keys(control), control)
    ]
    CallbackSet(callbacks...)
end


ame_trajectory_build_callback(
    control::Union{AMETrajectoryControl,AMETrajectoryDEControl},
) = build_callback(control)


ame_trajectory_build_callback(
    control::Union{AMETrajectoryControl,AMETrajectoryDEControl},
    sym::Symbol,
) = build_callback(control, sym)


ame_trajectory_build_callback(
    control::Union{InstPulseControl,InstDEPulseControl},
    sym::Symbol,
) = build_callback(control, sym, (c, pulse) -> c .= pulse * c)


ame_trajectory_build_callback(
    control::Union{FluctuatorControl,FluctuatorDEControl},
    sym::Symbol,
) = build_callback(control, sym)


function ame_trajectory_build_ode_function(H, ctr)
    update_func = build_update_func_ame_trajectory(ctr)
    cache = get_cache(H)
    diff_op = DiffEqArrayOperator(cache, update_func = update_func)
    jac_cache = similar(cache)
    jac_op = DiffEqArrayOperator(jac_cache, update_func = update_func)
    ff = ODEFunction(diff_op; jac_prototype = jac_op)
end


function build_update_func_ame_trajectory(
    ctr::Union{AMETrajectoryControl,AMETrajectoryDEControl},
)
    (A, u, p, t) -> update_cache!(A, p.control.op, p.tf, t)
end


function build_update_func_ame_trajectory(ctr::ControlSet)
    if has_fluctuator_control(ctr) && need_de_array(ctr, :fluctuator_control)
        update_func = function (A, u, p, t)
            update_cache!(A, p.control.ame_trajectory_control.op, p.tf, t)
            p.opensys(A, u, p.tf, t)
        end
    elseif has_fluctuator_control(ctr)
        update_func = function (A, u, p, t)
            update_cache!(A, p.control.ame_trajectory_control.op, p.tf, t)
            p.opensys(A, p.control.fluctuator_control(), p.tf, t)
        end
    else
        (A, u, p, t) ->
            update_cache!(A, p.control.ame_trajectory_control.op, p.tf, t)
    end
end
