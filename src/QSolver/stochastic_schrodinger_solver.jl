function build_ensemble_stochastic(
    A::Annealing,
    tf::Real,
    output_func,
    reduction;
    tspan=(0.0, tf),
    initializer=DEFAULT_INITIALIZER,
    kwargs...
)
    u0 = build_u0(A.u0, :v)
    initializer =
        initializer == DEFAULT_INITIALIZER ? (x, y) -> rand([-1, 1], x, y) :
        initializer

    flist = OpenQuantumBase.fluctuator_from_interactions(A.interactions)
    if length(flist) == 1
        cb = FluctuatorCallback(flist[1], initializer)
    else
        cb = CallbackSet([FluctuatorCallback(f, initializer) for f in flist]...)
    end
    ff = DiffEqLiouvillian(A.H, [], flist, size(A.H,1))
    p = ODEParams(ff, float(tf), A.annealing_parameter)

    update_func = function (cache, u, p, t)
        update_cache!(cache, p.L, p, t)
    end

    cache = get_cache(A.H)
    diff_op = DiffEqArrayOperator(cache, update_func=update_func)
    jac_cache = similar(cache)
    jac_op = DiffEqArrayOperator(jac_cache, update_func=update_func)
    ff = ODEFunction(diff_op, jac_prototype=jac_op)
    prob = ODEProblem{true}(ff, u0, tspan, p, callback=cb)

    ensemble_prob =
        EnsembleProblem(prob; output_func=output_func, reduction=reduction, kwargs...)
end
