"""
    function solve_redfield(
        A::Annealing,
        tf::Real,
        unitary;
        vectorize::Bool = false,
        dimensionless_time::Bool = true,
        tstops = Float64[],
        positivity_check::Bool = false,
        de_array_constructor = nothing,
        int_atol = 1e-8,
        int_rtol = 1e-6,
        kwargs...,
    )

Solve the time dependent Redfield equation for `Annealing` defined by `A` with total annealing time `tf`.

...
# Arguments
- `A::Annealing`: the Annealing object.
- `tf::Real`: the total annealing time.
- `unitary`: precalculated unitary of close system evolution.
- `vectorize::Bool = false`: whether to vectorize the density matrix.
- `dimensionless_time::Bool=true`: flag variable which, when set to true, informs the solver to work with dimensionless time.
- `tstops`: extra times that the timestepping algorithm must step to.
- `positivity_check::Bool = false`: whether to check the positivity of density matrix at each time step.
- `de_array_constructor = nothing`: the converting function if using `DEDataArray` type.
- `int_atol = 1e-8`: the absolute error tolerance for integration.
- `int_rtol = 1e-6`: the relative error tolerance for integration.
- `kwargs`: other keyword arguments supported by DifferentialEquations.jl.
...
"""
function solve_redfield(
    A::Annealing,
    tf::Real,
    unitary;
    vectorize::Bool = false,
    dimensionless_time::Bool = true,
    tstops = Float64[],
    positivity_check::Bool = false,
    de_array_constructor = nothing,
    int_atol = 1e-8,
    int_rtol = 1e-6,
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
    ff = __re_build_ode_function(A.H, vectorize)
    opensys = build_redfield(
        A.interactions,
        unitary,
        tf,
        atol = int_atol,
        rtol = int_rtol,
    )
    reset!(A.control)
    callback = __re_build_callback(A.control)
    p = ODEParams(A.H, tf; opensys = opensys, control = A.control)
    if positivity_check
        positivity_check_callback = FunctionCallingCallback(
            positivity_check_affect,
            func_everystep = true,
            func_start = false,
        )
        callback = CallbackSet(callback, positivity_check_callback)
    end
    prob = ODEProblem(ff, u0, (p) -> scaling_tspan(p.tf, A.sspan), p)
    solve(
        prob;
        alg_hints = [:nonstiff],
        tstops = tstops,
        callback = callback,
        kwargs...,
    )
end

"""
    function solve_CGME(
        A::Annealing,
        tf::Real,
        unitary;
        vectorize::Bool = false,
        dimensionless_time::Bool = true,
        tstops = Float64[],
        de_array_constructor = nothing,
        Ta = nothing,
        int_atol = 1e-8,
        int_rtol = 1e-6,
        kwargs...,
    )

Solve the time dependent CGME for `Annealing` defined by `A` with total annealing time `tf`.

...
# Arguments
- `A::Annealing`: the Annealing object.
- `tf::Real`: the total annealing time.
- `unitary`: precalculated unitary of close system evolution.
- `vectorize::Bool = false`: whether to vectorize the density matrix.
- `dimensionless_time::Bool=true`: flag variable which, when set to true, informs the solver to work with dimensionless time.
- `tstops`: extra times that the timestepping algorithm must step to.
- `de_array_constructor = nothing`: the converting function if using `DEDataArray` type.
- `Ta = nothing`: coarse-graining time.
- `int_atol = 1e-8`: the absolute error tolerance for integration.
- `int_rtol = 1e-6`: the relative error tolerance for integration.
- `kwargs`: other keyword arguments supported by DifferentialEquations.jl.
...
"""
function solve_CGME(
    A::Annealing,
    tf::Real,
    unitary;
    vectorize::Bool = false,
    dimensionless_time::Bool = true,
    tstops = Float64[],
    de_array_constructor = nothing,
    Ta = nothing,
    int_atol = 1e-8,
    int_rtol = 1e-6,
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
    ff = __re_build_ode_function(A.H, vectorize)
    opensys = build_CGME(
        A.interactions,
        unitary,
        tf,
        atol = int_atol,
        rtol = int_rtol,
        Ta = Ta,
    )
    reset!(A.control)
    callback = __re_build_callback(A.control)
    p = ODEParams(A.H, tf; opensys = opensys, control = A.control)
    prob = ODEProblem(ff, u0, (p) -> scaling_tspan(p.tf, A.sspan), p)
    solve(
        prob;
        alg_hints = [:nonstiff],
        tstops = tstops,
        callback = callback,
        kwargs...,
    )
end

__re_build_callback(control) = nothing
__re_build_callback(control::Union{InstPulseControl,InstDEPulseControl}) =
    build_callback(control, pulse_on_density!)

function __re_build_ode_function(H, vectorize)
    if vectorize == false
        redfield_f
    else
        cache = get_cache(H, vectorize)
        diff_op =
            DiffEqArrayOperator(cache, update_func = redfield_vectorize_update!)
        jac_cache = similar(cache)
        jac_op = DiffEqArrayOperator(
            jac_cache,
            update_func = redfield_vectorize_update!,
        )
        ff = ODEFunction(diff_op; jac_prototype = jac_op)
    end
end


function redfield_f(du, u, p, t)
    p.H(du, u, p.tf, t)
    p.opensys(du, u, p.tf, t)
end


function redfield_vectorize_update!(A, u, p, t)
    update_vectorized_cache!(A, p.H, p.tf, t)
    update_vectorized_cache!(A, p.opensys, p.tf, t)
end


# ================ the following codes are for hybrid Redfield =================
function build_ensemble_hybrid_redfield(
    A::Annealing,
    tf::Real,
    unitary,
    output_func,
    reduction;
    vectorize::Bool = false,
    dimensionless_time::Bool = true,
    positivity_check::Bool = false,
    de_array_constructor = nothing,
    fluctuator_de_field = nothing,
    int_atol = 1e-8,
    int_rtol = 1e-6,
    initializer = DEFAULT_INITIALIZER,
    tstops = Float64[],
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
        needed_symbol = fluctuator_de_field == nothing ? [] :
                            [fluctuator_de_field],
    )

    control, fluctuator_opensys, redfield_opensys =
        build_hybrid_redfield_control_from_interactions(
            A.interactions,
            unitary,
            tf,
            int_atol,
            int_rtol,
            fluctuator_de_field,
        )
    control = A.control == nothing ? control : ControlSet(control, A.control)
    p = ODEParams(tf; control = control, opensys = opensys)
    callback = __rehybrid_build_callback(control)
    if positivity_check
        positivity_check_callback = FunctionCallingCallback(
            positivity_check_affect,
            func_everystep = true,
            func_start = false,
        )
        callback = CallbackSet(callback, positivity_check_callback)
    end

    ff = __rehybrid_build_ode_function(A.H, control)
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

function __rehybrid_build_callback(control::ControlSet)
    callbacks = [
        __rehybrid_build_callback(v, k)
        for (k, v) in zip(keys(control), control)
    ]
    CallbackSet(callbacks...)
end

__rehybrid_build_callback(
    control::Union{InstPulseControl,InstDEPulseControl},
    sym::Symbol,
) = build_callback(control, sym, pulse_on_density!)

__rehybrid_build_callback(
    control::Union{FluctuatorControl,FluctuatorDEControl},
    sym::Symbol,
) = build_callback(control, sym)
