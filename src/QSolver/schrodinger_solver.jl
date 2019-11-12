function solve_schrodinger(A::Annealing, tf::Real; span_unit = false, tstops = Float64[], kwargs...)
    tf = prepare_tf(tf, span_unit)
    tstops = prepare_tstops(tf, tstops, A.tstops)
    u0 = prepare_u0(A.u0, type = :v, control = A.control)
    p = AnnealingParams(A.H, tf; control = A.control)
    callback = construct_callback(A.control, :schrodinger)
    ff = schrodinger_construct_ode_function(A.H, A.control)
    prob = ODEProblem{true}(ff, u0, (p) -> scaling_tspan(p.tf, A.sspan), p)
    solve(prob; alg_hints = [:nonstiff], callback = callback, tstops = tstops, kwargs...)
end


function solve_schrodinger(
    A::Annealing,
    tf::Vector{T},
    alg,
    para_alg = EnsembleSerial();
    output_func = (sol, i) -> (sol, false),
    span_unit = false,
    tstops = Float64[],
    kwargs...,
) where {T<:Real}
    u0 = prepare_u0(A.u0, type = :v, control = A.control)
    # the hard code 1.0 is just to set the type of argument correct
    t0 = prepare_tf(1.0, span_unit)
    tstops = prepare_tstops(1.0, tstops, A.tstops)
    p = AnnealingParams(A.H, t0; control = A.control)
    # set the type of tf array
    tf_arr = float.(tf)
    # resolve control
    callback = construct_callback(A.control, :schrodinger)
    ff = schrodinger_construct_ode_function(A.H, A.control)
    #
    prob_func = (prob, i, repeat) -> begin
        p = set_tf(prob.p, tf_arr[i])
        ODEProblem{true}(prob.f, prob.u0, prob.tspan, p)
    end
    if span_unit == true
        if !isempty(tstops)
            callback = CallbackSet(callback, PresetTimeCallback(tstops, empty_affect!, save_positions=(true, false)))
        end
        tstops = []
    end
    prob = ODEProblem{true}(ff, u0, (p) -> scaling_tspan(p.tf, A.sspan), p)
    ensemble_prob = EnsembleProblem(prob; prob_func = prob_func, output_func = output_func)
    solve(
        ensemble_prob,
        alg,
        para_alg;
        trajectories = length(tf),
        tstops = tstops,
        callback = callback,
        kwargs...,
    )
end


function schrodinger_construct_ode_function(H, ::Union{Nothing,InstPulseControl})
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
