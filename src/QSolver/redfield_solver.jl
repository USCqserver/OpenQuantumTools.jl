"""
$(SIGNATURES)

Solve the time dependent Redfield equation for `Annealing` defined by `A` with total annealing time `tf`.

...
# Arguments
- `A::Annealing`: the Annealing object.
- `tf::Real`: the total annealing time.
- `unitary`: precalculated unitary of close system evolution.
- `vectorize::Bool = false`: whether to vectorize the density matrix.
- `int_atol = 1e-8`: the absolute error tolerance for integration.
- `int_rtol = 1e-6`: the relative error tolerance for integration.
- `Ta = tf`: the time scale for backward integration.
- `kwargs`: other keyword arguments supported by DifferentialEquations.jl.
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
    L = build_redfield(A.interactions, unitary, Ta, int_atol, int_rtol)
    R = RedfieldOperator(A.H, L)

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

Solve the time dependent CGME for `Annealing` defined by `A` with total annealing time `tf`.

...
# Arguments
- `A::Annealing`: the Annealing object.
- `tf::Real`: the total annealing time.
- `unitary`: precalculated unitary of close system evolution.
- `vectorize::Bool = false`: whether to vectorize the density matrix.
- `Ta = nothing`: coarse-graining time. If set to nothing, the solver will automatically choose the value.
- `int_atol = 1e-8`: the absolute error tolerance for integration.
- `int_rtol = 1e-6`: the relative error tolerance for integration.
- `kwargs`: other keyword arguments supported by DifferentialEquations.jl.
...
"""
function solve_CGME(
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
    L = build_CGG(A.interactions, unitary, tf, Ta,
            int_atol, int_rtol)
    R = RedfieldOperator(A.H, L)
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
    reds, stocs = build_red_lvs(A.interactions, unitary, Ta, int_atol, int_rtol)
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
    R = RedfieldOperator(A.H, [reds; stocs])
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

function build_red_lvs(iset, U, Ta, atol, rtol)
    reds = []
    stochastic = []
    for i in iset
        if typeof(i.bath) <: EnsembleFluctuator
            push!(stochastic, build_fluctuator(i.coupling, i.bath))
        else
            push!(reds, build_redfield(i.coupling, i.bath, U, Ta, atol, rtol))
        end
    end
    reds, stochastic
end
