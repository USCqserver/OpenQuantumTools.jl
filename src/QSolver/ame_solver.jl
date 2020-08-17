"""
$(SIGNATURES)

Solve the adiabatic master equation for `Annealing` defined by `A` with total annealing time `tf`.

...
# Arguments
- `A::Annealing`: the Annealing object.
- `tf::Real`: the total annealing time.
- `tspan` = (0, tf): time interval to solve.
- `ω_hint=[]` : grid for precalculating the lambshift; skip the precalculation if empty.
- 'lambshift::Bool=true' : whether to include Lambshift in the calculation.
- `lvl::Int=size(A.H, 1)` : number of levels to keep. The default value is the dimension for the Hamiltonian.
- `vectorize::Bool = false`: whether to vectorize the density matrix.
- `kwargs` : other keyword arguments supported by DifferentialEquations.jl.
...
"""
function solve_ame(
    A::Annealing,
    tf::Real;
    tspan=(0.0, tf),
    ω_hint=[],
    lambshift::Bool=true,
    lvl::Int=size(A.H, 1),
    vectorize::Bool=false,
    kwargs...,
)
    u0 = build_u0(A.u0, :m, vectorize=vectorize)
    davies = build_davies(A.interactions, ω_hint, lambshift)
    if vectorize
        error("Vectorization is not yet supported for adiabatic master equation.")
    end
    f = AMEOperator(A.H, davies, lvl)
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
    lvl::Int=size(A.H, 1),
    initializer=initializer,
    kwargs...,
)
    u0 = build_u0(A.u0, :v)

    dlist, flist = build_ametr_lvs(A.interactions, ω_hint, lambshift)
    L, cb = build_ametr_op_callback(A.H, lvl, dlist, flist, initializer)
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

function build_ametr_lvs(iset, ω_hint, lambshift)
    davies = []
    stochastic = []
    for i in iset
        if typeof(i.bath) <: EnsembleFluctuator
            push!(stochastic, build_fluctuator(i.coupling, i.bath))
        else
            push!(davies, build_davies(i.coupling, i.bath, ω_hint, lambshift))
        end
    end
    davies, stochastic
end

function build_ametr_op_callback(H, lvl, dlist, flist, initializer)
    initializer =
        initializer == DEFAULT_INITIALIZER ? (x, y) -> rand([-1, 1], x, y) :
        initializer
    if isempty(flist)
        f = AMEOperator(H, dlist, lvl)
        cb = AMEtrajectoryCallback()
    elseif isempty(dlist)
        error("No interactions support AME. Use other ensemble instead.")
    else
        f = QTBase.OpenSysOpHybrid(H, dlist, flist, lvl)
        cb = CallbackSet(
            AMEtrajectoryCallback(),
            [FluctuatorCallback(f, initializer) for f in flist]...,
        )
    end
    f, cb
end
