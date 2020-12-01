"""
$(SIGNATURES)

Solve the Lindblad equation defined by `A` for a total evolution time `tf`.

...
# Arguments
- `A::Annealing`: the Annealing object.
- `tf::Real`: the total annealing time.
- `tspan` = (0, tf): time interval to solve the dynamics.
- `vectorize::Bool = false`: whether to vectorize the density matrix.
- `kwargs` : other keyword arguments supported by `DifferentialEquations`.
...
"""
function solve_lindblad(
    A::Annealing,
    tf::Real;
    tspan=(0.0, tf),
    vectorize::Bool=false,
    kwargs...,
)
    u0 = build_u0(A.u0, :m, vectorize=vectorize)
    L = OpenQuantumBase.lindblad_from_interactions(A.interactions)
    if vectorize
        error("Vectorization is not yet supported for Lindblad equation.")
    end
    f = DiffEqLiouvillian(A.H, [], L, size(A.H, 1))
    p = ODEParams(f, float(tf), A.annealing_parameter)
    prob = ODEProblem(f, u0, tspan, p)
    solve(prob; alg_hints=[:nonstiff], kwargs...)
end

function build_ensemble_lindblad(
    A::Annealing,
    tf::Real,
    output_func,
    reduction;
    tspan=(0.0, tf),
    kwargs...,
)
    u0 = build_u0(A.u0, :v)
    L = DiffEqLiouvillian(A.H, [], OpenQuantumBase.lindblad_from_interactions(A.interactions), size(A.H, 1))
    cb = LindbladJumpCallback()
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