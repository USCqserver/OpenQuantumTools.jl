"""
$(SIGNATURES)

Solve the time-dependent Redfield equation defined by `A` for a total evolution time `tf`.

...
# Arguments
- `A::Annealing`: the `Annealing`/`Evolution` object.
- `tf::Real`: the total evolution time.
- `unitary`: precomputed unitary operator of the corresponding closed-system evolution.
- `vectorize::Bool = false`: whether to vectorize the density matrix.
- `int_atol = 1e-8`: the absolute error tolerance for integration.
- `int_rtol = 1e-6`: the relative error tolerance for integration.
- `Ta = tf`: the timescale of the backward integration.
- `kwargs`: other keyword arguments supported by `DifferentialEquations`.
...
"""
function solve_redfield(
    A::Annealing,
    tf::Real,
    unitary;
    vectorize::Bool=false,
    int_atol=1e-8,
    int_rtol=1e-6,
    Ta=tf,
    kwargs...,
)
    u0 = build_u0(A.u0, :m, vectorize=vectorize)
    L = OpenQuantumBase.redfield_from_interactions(A.interactions, unitary, Ta, int_atol, int_rtol)
    R = DiffEqLiouvillian(A.H, [], L, size(A.H, 1))

    update_func = function (A, u, p, t)
        update_vectorized_cache!(A, p.L, p, t)
    end

    if vectorize == false
        ff = ODEFunction{true}(R)
    else
        cache = vectorize_cache(get_cache(A.H))
        diff_op = DiffEqArrayOperator(cache, update_func=update_func)
        jac_cache = similar(cache)
        jac_op = DiffEqArrayOperator(jac_cache, update_func=update_func)
        ff = ODEFunction(diff_op, jac_prototype=jac_op)
    end

    p = ODEParams(R, float(tf), A.annealing_parameter)
    prob = ODEProblem(ff, u0, (0.0, float(tf)), p)
    solve(prob; alg_hints=[:nonstiff], kwargs...)
end

"""
$(SIGNATURES)

Solve the time-dependent coarse-grained master equation (CGME) defined by `A` for a total evolution time `tf`.

...
# Arguments
- `A::Annealing`: the `Annealing`/`Evolution` object.
- `tf::Real`: the total evolution time.
- `unitary`: precomputed unitary operator of the corresponding closed-system evolution.
- `vectorize::Bool = false`: whether to vectorize the density matrix.
- `Ta = nothing`: coarse-graining time. If set to nothing, the solver will automatically choose the value.
- `int_atol = 1e-8`: the absolute error tolerance for integration.
- `int_rtol = 1e-6`: the relative error tolerance for integration.
- `kwargs`: other keyword arguments supported by `DifferentialEquations`.
...
"""
function solve_cgme(
    A::Annealing,
    tf::Real,
    unitary;
    vectorize::Bool=false,
    Ta=nothing,
    int_atol=1e-8,
    int_rtol=1e-6,
    kwargs...,
)
    u0 = build_u0(A.u0, :m, vectorize=vectorize)
    L = OpenQuantumBase.cg_from_interactions(A.interactions, unitary, tf, Ta,
            int_atol, int_rtol)
    R = DiffEqLiouvillian(A.H, [], L, size(A.H, 1))
    update_func = function (A, u, p, t)
        update_vectorized_cache!(A, p.L, p, t)
    end

    if vectorize == false
        ff = ODEFunction{true}(R)
    else
        cache = vectorize_cache(get_cache(A.H))
        diff_op = DiffEqArrayOperator(cache, update_func=update_func)
        jac_cache = similar(cache)
        jac_op = DiffEqArrayOperator(jac_cache, update_func=update_func)
        ff = ODEFunction(diff_op, jac_prototype=jac_op)
    end

    p = ODEParams(R, float(tf), A.annealing_parameter)
    prob = ODEProblem(ff, u0, (0.0, float(tf)), p)
    solve(prob; alg_hints=[:nonstiff], kwargs...)
end

"""
$(SIGNATURES)

Solve the time-dependent universal Lindblad equation (ULE) defined by `A` for a total evolution time `tf`.

...
# Arguments
- `A::Annealing`: the `Annealing`/`Evolution` object.
- `tf::Real`: the total evolution time.
- `unitary`: precomputed unitary operator of the corresponding closed-system evolution.
- `vectorize::Bool = false`: whether to vectorize the density matrix.
- `int_atol = 1e-8`: the absolute error tolerance for integration.
- `int_rtol = 1e-6`: the relative error tolerance for integration.
- `Ta = tf`: the timescale of the integration region.
- `kwargs`: other keyword arguments supported by `DifferentialEquations`.
...
"""
function solve_ule(
    A::Annealing,
    tf::Real,
    unitary;
    vectorize::Bool=false,
    int_atol=1e-8,
    int_rtol=1e-6,
    Ta=tf,
    kwargs...,
)
    u0 = build_u0(A.u0, :m, vectorize=vectorize)
    L = OpenQuantumBase.ule_from_interactions(A.interactions, unitary, Ta, int_atol, int_rtol)
    R = DiffEqLiouvillian(A.H, [], L, size(A.H, 1))

    update_func = function (A, u, p, t)
        update_vectorized_cache!(A, p.L, p, t)
        error("Vectorization of ULE is not currently supported.")
    end

    if vectorize == false
        ff = ODEFunction{true}(R)
    else
        cache = vectorize_cache(get_cache(A.H))
        diff_op = DiffEqArrayOperator(cache, update_func=update_func)
        jac_cache = similar(cache)
        jac_op = DiffEqArrayOperator(jac_cache, update_func=update_func)
        ff = ODEFunction(diff_op, jac_prototype=jac_op)
    end

    p = ODEParams(R, float(tf), A.annealing_parameter)
    prob = ODEProblem(ff, u0, (0.0, float(tf)), p)
    solve(prob; alg_hints=[:nonstiff], kwargs...)
end

function build_ensemble_redfield(
    A::Annealing,
    tf,
    unitary,
    output_func,
    reduction;
    vectorize::Bool=false,
    int_atol=1e-8,
    int_rtol=1e-6,
    Ta=tf,
    initializer=DEFAULT_INITIALIZER,
    kwargs...,
)
    u0 = build_u0(A.u0, :m, vectorize=vectorize)
    stocs = OpenQuantumBase.fluctuator_from_interactions(A.interactions)
    reds = OpenQuantumBase.redfield_from_interactions(A.interactions, unitary, Ta, int_atol, int_rtol)
    if isempty(stocs)
        error("No stochastic bath detected. Use the normal Redfeidl solver instead.")
    end
    # set the initializer
    initializer =
        initializer == DEFAULT_INITIALIZER ? (x, y) -> rand([-1, 1], x, y) :
        initializer

    if length(stocs) == 1
        cb = FluctuatorCallback(stocs[1], initializer)
    else
        cb = CallbackSet([FluctuatorCallback(f, initializer) for f in stocs]...)
    end
    R = DiffEqLiouvillian(A.H, [], [reds; stocs], size(A.H, 1))
    p = ODEParams(R, float(tf), A.annealing_parameter)

    update_func = function (A, u, p, t)
        update_vectorized_cache!(A, p.L, p, t)
    end

    if vectorize == false
        ff = ODEFunction{true}(R)
    else
        cache = vectorize_cache(get_cache(A.H))
        diff_op = DiffEqArrayOperator(cache, update_func=update_func)
        jac_cache = similar(cache)
        jac_op = DiffEqArrayOperator(jac_cache, update_func=update_func)
        ff = ODEFunction(diff_op, jac_prototype=jac_op)
    end

    prob = ODEProblem{true}(ff, u0, (0.0, float(tf)), p, callback=cb)

    ensemble_prob =
        EnsembleProblem(prob; output_func=output_func, reduction=reduction)

    ensemble_prob
end
