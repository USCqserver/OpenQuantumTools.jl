"""
$(SIGNATURES)

Solve the adiabatic master equation defined by `A` for a total evolution time `tf`.

...
# Arguments
- `A::Annealing`: the `Annealing`/`Evolution` object.
- `tf::Real`: the total annealing time.
- `tspan = (0, tf)`: time interval to solve the dynamics.
- `ω_hint=[]`: specify a grid to precompute the ``S`` function in the Lamb shift term. Skip the precomputation if empty.
- `digits::Int=8: the number of digits to keep when checking if a gap is zero.`
- `sigdigits::Int=8: the number of significant digits when rounding non-zero gaps for comparison.`
- `lambshift::Bool=true`: whether to include the Lamb shift in the simulation.
- `lambshift_kwargs=Dict()`: keyword arguments to be supplied to the custom routine to calculate the `S` function in the Lamb shift term.
- `lvl::Int=size(A.H, 1)`: number of levels to keep. The higher levels are ignored to speed up the computation.
- `vectorize::Bool = false`: whether to vectorize the density matrix.
- `one_sided::Bool=false`: whether to solve the one-sided AME.
- `kwargs` : other keyword arguments supported by `DifferentialEquations`.
...
"""
function solve_ame(
    A::Annealing{false},
    tf::Real;
    tspan=(0.0, tf),
    ω_hint=[],
    lambshift::Bool=true,
    lambshift_kwargs=Dict(),
    lvl::Int=size(A.H, 1),
    vectorize::Bool=false,
    one_sided::Bool=false,
    digits::Int=8,
    sigdigits::Int=8,
    kwargs...
)
    u0 = build_u0(A.u0, :m, vectorize=vectorize)
    if one_sided == false
        L = OpenQuantumBase.davies_from_interactions(A.interactions, ω_hint, lambshift, lambshift_kwargs)
    else
        L = OpenQuantumBase.onesided_ame_from_interactions(A.interactions, ω_hint, lambshift, lambshift_kwargs)
    end
    if vectorize
        throw(ArgumentError("Vectorization is not yet supported for adiabatic master equation."))
    end
    f = DiffEqLiouvillian(A.H, L, [], lvl, digits=digits, sigdigits=sigdigits)
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
    lambshift_kwargs=nothing,
    lvl::Int=size(A.H, 1),
    initializer=initializer,
    save_positions=(false, false),
    digits::Int=8,
    sigdigits::Int=8,
    kwargs...
)
    u0 = build_u0(A.u0, :v)

    flist = OpenQuantumBase.fluctuator_from_interactions(A.interactions)
    dlist = OpenQuantumBase.davies_from_interactions(A.interactions, ω_hint, lambshift, lambshift_kwargs)
    L = DiffEqLiouvillian(A.H, dlist, flist, lvl, digits=digits, sigdigits=sigdigits)
    cb = build_jump_callback(dlist, flist, initializer, save_positions)
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

function build_jump_callback(dlist, flist, initializer, save_positions)
    initializer =
        initializer == DEFAULT_INITIALIZER ? (x, y) -> rand([-1, 1], x, y) :
        initializer
    if isempty(flist)
        cb = LindbladJumpCallback(save_positions)
    elseif isempty(dlist)
        error("No interactions support AME. Use other ensemble instead.")
    else
        cb = CallbackSet(
            LindbladJumpCallback(save_positions),
            [FluctuatorCallback(f, initializer, save_positions) for f in flist]...,
        )
    end
    cb
end

function solve_ame(
    A::Annealing{true},
    tf::Real;
    tspan=(0.0, tf),
    ω_hint=[],
    lambshift::Bool=true,
    lambshift_kwargs=Dict(),
    lvl::Int=size(A.H, 1),
    vectorize::Bool=false,
    one_sided::Bool=false,
    digits::Int=8,
    sigdigits::Int=8,
    cutoff::Real=Inf,
    kwargs...
)
    # rotate the system into the eigen state of the Hamiltonian
    w, v = eigen_decomp(A.H, lvl=lvl)
    H = w |> Diagonal |> sparse |> Hamiltonian
    inters = rotate(A.interactions, v)
    # build gaps for AME
    gap_idx = OpenQuantumBase.build_gap_indices(w, digits, sigdigits, cutoff, lvl)
    # prepare for the initial state
    u0 = build_u0(A.u0, :m, vectorize=vectorize)
    u0 = v' * u0 * v

    if one_sided == false
        L = OpenQuantumBase.davies_from_interactions(gap_idx, inters, ω_hint, lambshift, lambshift_kwargs)
    else
        throw(ArgumentError("One-sided AME not supported for constant Hamiltonian type."))
    end
    if vectorize
        throw(ArgumentError("Vectorization is not yet supported for adiabatic master equation."))
    end
    f = OpenQuantumBase.build_diffeq_liouvillian(H, [], L, lvl, digits=digits, sigdigits=sigdigits)
    p = ODEParams(f, float(tf), A.annealing_parameter)
    prob = ODEProblem(f, u0, tspan, p)
    solve(prob; alg_hints=[:nonstiff], kwargs...)
end
