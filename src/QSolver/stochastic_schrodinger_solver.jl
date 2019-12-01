function solve_stochastic_schrodinger(
    A::Annealing,
    tf::Real,
    trajectories::Int,
    alg,
    para_alg = EnsembleSerial();
    output_func = (sol, i) -> (sol, false),
    span_unit = false,
    tstops = Float64[],
    kwargs...,
)
    ensemble_prob, kwarg_dict = build_ensemble_problem(
        A,
        tf,
        :stochastic_schrodinger,
        output_func = output_func,
        span_unit = span_unit,
        tstops = tstops;
        kwargs...,
    )
    solve(
        ensemble_prob,
        alg,
        para_alg;
        trajectories = trajectories,
        kwarg_dict...,
    )
end


function build_ensemble_problem_stochastic_schrodinger(
    A::Annealing,
    tf::Real;
    output_func = (sol, i) -> (sol, false),
    span_unit = false,
    tstops = Float64[],
    kwargs...,
)
    tf = build_tf(tf, span_unit)
    tstops = build_tstops(tf, tstops, A.tstops)
    # build control object from bath; a prototype implementation
    control = FluctuatorControl(tf, length(A.coupling), A.bath)
    control = A.control == nothing ? control : ControlSet(control, A.control)
    u0 = build_u0(A.u0, :v, control = control)
    opensys = build_stochastic_opensys(A.coupling)
    callback = build_callback(control, :stochastic_schrodinger)
    p = AnnealingParams(A.H, tf; opensys = opensys, control = control)
    ff = stochastic_schrodinger_build_ode_function(A.H, A.control)
    prob = ODEProblem{true}(ff, u0, (p) -> scaling_tspan(p.tf, A.sspan), p)

    ensemble_prob = EnsembleProblem(
        prob;
        prob_func = DEFAULT_FLUCTUATOR_CONTROL_PROB_FUNC,
        output_func = output_func,
    )

    kwarg_dict = Dict(kwargs...)
    kwarg_dict[:tstops] = tstops
    kwarg_dict[:callback] = callback
    ensemble_prob, kwarg_dict
end


function stochastic_schrodinger_build_ode_function(H, control)
    cache = get_cache(H)
    diff_op = DiffEqArrayOperator(
        cache,
        update_func = (A, u, p, t) -> stochastic_update!(
            A,
            u,
            p.tf,
            t,
            p.H,
            p.opensys,
        ),
    )
    jac_cache = similar(cache)
    jac_op = DiffEqArrayOperator(
        jac_cache,
        update_func = (A, u, p, t) -> stochastic_update!(
            A,
            u,
            p.tf,
            t,
            p.H,
            p.opensys,
        ),
    )
    ff = ODEFunction(diff_op; jac_prototype = jac_op)
end


function stochastic_update!(A, u, tf, t, H, opensys)
    update_cache!(A, H, tf, t)
    opensys(A, u, tf, t)
end
