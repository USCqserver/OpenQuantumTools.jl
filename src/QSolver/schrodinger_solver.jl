function solve_schrodinger(
    A::Annealing,
    tf::Real;
    span_unit = false,
    tstops = Float64[],
    kwargs...,
)
    tf = build_tf(tf, span_unit)
    tstops = build_tstops(tf, tstops, A.tstops)
    u0 = build_u0(A.u0, :v, control = A.control)
    p = AnnealingParams(A.H, tf; control = A.control)
    callback = build_callback(A.control, :schrodinger)
    ff = schrodinger_construct_ode_function(A.H, A.control)
    prob = ODEProblem{true}(ff, u0, (p) -> scaling_tspan(p.tf, A.sspan), p)
    solve(
        prob;
        alg_hints = [:nonstiff],
        callback = callback,
        tstops = tstops,
        kwargs...,
    )
end


function build_ensemble_problem_schrodinger(
    A::Annealing,
    prob_func,
    output_func,
    reduction;
    span_unit = false,
) where T<:Number
    u0 = build_u0(A.u0, :v, control = A.control)
    # the hard code 1.0 is just to set the type of argument correct
    t0 = build_tf(1.0, span_unit)
    p = AnnealingParams(A.H, t0; control = A.control)
    # resolve control
    callback = build_callback(A.control, :schrodinger)
    ff = schrodinger_construct_ode_function(A.H, A.control)
    prob = ODEProblem{true}(ff, u0, (p) -> scaling_tspan(p.tf, A.sspan), p)
    ensemble_prob = EnsembleProblem(
        prob;
        prob_func = prob_func,
        output_func = output_func,
    )
    ensemble_prob, callback
end


function schrodinger_construct_ode_function(
    H,
    ::Union{Nothing,InstPulseControl},
)
    cache = get_cache(H)
    diff_op = DiffEqArrayOperator(
        cache,
        update_func = (A, u, p, t) -> update_cache!(A, p.H, p.tf, t),
    )
    jac_cache = similar(cache)
    jac_op = DiffEqArrayOperator(
        jac_cache,
        update_func = (A, u, p, t) -> update_cache!(A, p.H, p.tf, t),
    )
    ff = ODEFunction(diff_op; jac_prototype = jac_op)
end


function schrodinger_construct_ode_function(
    H,
    pause_control::PausingControl,
)
    cache = get_cache(H)
    diff_op = DiffEqArrayOperator(
        cache,
        update_func = (A, u, p, t) -> update_cache!(A, u, p.tf, t, p.H, p.control),
    )
    jac_cache = similar(cache)
    jac_op = DiffEqArrayOperator(
        jac_cache,
        update_func = (A, u, p, t) -> update_cache!(A, u, p.tf, t, p.H, p.control),
    )
    ff = ODEFunction(diff_op; jac_prototype = jac_op)
end
