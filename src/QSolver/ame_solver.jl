"""
$(SIGNATURES)

Solve the adiabatic master equation defined by `A` for a total evolution time `tf`.

...
# Arguments
- `A::Annealing`: the `Annealing`/`Evolution` object.
- `tf::Real`: the total annealing time.
- `tspan = (0, tf)`: time interval to solve the dynamics.
- `ω_hint=[]`: specify a grid to precompute the ``S`` function in the Lamb shift term. Skip the precomputation if empty.
- `lambshift::Bool=true`: whether to include the Lamb shift in the simulation.
- `lambshift_S=nothing`: provide a custom routine to calculate the `S` function in the Lamb shift term. 
- `lvl::Int=size(A.H, 1)`: number of levels to keep. The higher levels are ignored to speed up the computation.
- `vectorize::Bool = false`: whether to vectorize the density matrix.
- `one_sided=false`: whether to solve the one-sided AME.
- `kwargs` : other keyword arguments supported by `DifferentialEquations`.
...
"""
function solve_ame(
    A::Annealing,
    tf::Real;
    tspan=(0.0, tf),
    ω_hint=[],
    lambshift::Bool=true,
    lambshift_S=nothing,
    lvl::Int=size(A.H, 1),
    vectorize::Bool=false,
    one_sided=false,
    kwargs...,
)
    u0 = build_u0(A.u0, :m, vectorize=vectorize)
    if one_sided == false
        L = QTBase.davies_from_interactions(A.interactions, ω_hint, lambshift, lambshift_S)
    else
        L = build_onesided_ame(A.interactions, ω_hint, lambshift)
    end
    if vectorize
        error("Vectorization is not yet supported for adiabatic master equation.")
    end
    f = DiffEqLiouvillian(A.H, L, [], lvl)
    p = ODEParams(f, float(tf), A.annealing_parameter)
    prob = ODEProblem(f, u0, tspan, p)
    solve(prob; alg_hints=[:nonstiff], kwargs...)
end

function build_ensemble_ame(
    A::Annealing,
    tf::Real,
    output_func,
    reduction;
    tspan=(0.0, tf),
    ω_hint=[],
    lambshift::Bool=true,
    lambshift_S=nothing,
    lvl::Int=size(A.H, 1),
    initializer=initializer,
    kwargs...,
)
    u0 = build_u0(A.u0, :v)

    flist = QTBase.fluctuator_from_interactions(A.interactions)
    dlist = QTBase.davies_from_interactions(A.interactions, ω_hint, lambshift, lambshift_S)
    L = DiffEqLiouvillian(A.H, dlist, flist, lvl)
    cb = build_jump_callback(dlist, flist, initializer)
    p = ODEParams(L, tf, A.annealing_parameter)
    update_func = function (cache, u, p, t)
        update_cache!(cache, p.L, p, t)
    end
    cache = get_cache(A.H)
    diff_op = DiffEqArrayOperator(cache, update_func=update_func)
    jac_cache = similar(cache)
    jac_op = DiffEqArrayOperator(jac_cache, update_func=update_func)
    ff = ODEFunction(diff_op; jac_prototype=jac_op)

    prob = ODEProblem{true}(ff, u0, tspan, p, callback=cb)
    ensemble_prob = EnsembleProblem(
        prob;
        output_func=output_func,
        reduction=reduction,
        kwargs...
    )
    ensemble_prob
end

function build_jump_callback(dlist, flist, initializer)
    initializer =
        initializer == DEFAULT_INITIALIZER ? (x, y) -> rand([-1, 1], x, y) :
        initializer
    if isempty(flist)
        cb = LindbladJumpCallback()
    elseif isempty(dlist)
        error("No interactions support AME. Use other ensemble instead.")
    else
        cb = CallbackSet(
            LindbladJumpCallback(),
            [FluctuatorCallback(f, initializer) for f in flist]...,
        )
    end
    cb
end
