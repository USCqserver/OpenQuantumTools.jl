"""
    function solve_ame(
        A::Annealing,
        tf::Real;
        span_unit::Bool = false,
        ω_hint = [],
        lambshift::Bool = true,
        lvl,
        eig_tol = 1e-8,
        maxiter = 3000,
        ncv,
        v0,
        kwargs...,
    )

Solve the adiabatic master equation for `Annealing` defined by `A` with total annealing time `tf`.

...
# Arguments
- `A::Annealing`: the Annealing object.
- `tf::Real`: the total annealing time.
- `span_unit::Bool=false`: flag variable which, when set to true, informs the solver to work with time in physical unit.
- `ω_hint=[]` : grid for precalculating the lambshift; skip the precalculation if empty.
- 'lambshift::Bool=true' : whether to include Lambshift in the calculation.
- `lvl::Int` : number of levels to keep. The default value is the dim/dim-1 for dense/sparse Hamiltonian.
- `eig_tol = 1e-8` : (only for sparse Hamiltonian) relative tolerance for eigen-decomposition.
- 'ncv' : (only for sparse Hamiltonian) number of Krylov vectors used in Arpack. The default value is max(20, 2*lvl+1).
- `maxiter::Int = 3000` : (only for sparse Hamiltonian) maximum number of iterations for Arpack.
- `v0` : (only for sparse Hamiltonian) starting vector for Arpack. The default value is the initial state `u0` if it is a state vector and ground state of the initial Hamiltonian if not.
- `kwargs` : other keyword arguments supported by DifferentialEquations.jl
...
"""
function solve_ame(
    A::Annealing,
    tf::Real;
    span_unit::Bool = false,
    ω_hint = [],
    lambshift::Bool = true,
    lvl::Int = size(A.H, 1),
    eig_tol = 0.0,
    maxiter::Int = 3000,
    ncv = 20,
    v0 = copy(A.u0),
    kwargs...,
)
    u0 = build_u0(A.u0, :m, control = A.control)
    tf = build_tf(tf, span_unit)
    davies = build_davies(A.coupling, A.bath, ω_hint, lambshift)
    f = AMEDiffEqOperator(
        A.H,
        davies,
        lvl,
        eig_tol = eig_tol,
        control = A.control,
        maxiter = maxiter,
        ncv = ncv,
        v0 = v0,
    )
    p = LightAnnealingParams(tf; control = A.control)
    if typeof(A.control) <: PausingControl
        cb = DiscreteCallback(pause_condition, pause_affect!)
        kwargs = Dict{Symbol,Any}(kwargs)
        kwargs[:callback] = cb
    end
    tspan, tstops = scaling_time(tf, A.sspan, A.tstops)
    prob = ODEProblem(f, u0, tspan, p)
    solve(prob; alg_hints = [:nonstiff], tstops = tstops, kwargs...)
end


function (D::QTBase.AMEDenseDiffEqOperator{true,T})(du, u, p, t) where {T<:PausingControl}
    s, a_scale, g_scale = p.control(p.tf, t)
    hmat = D.H(u, a_scale, g_scale, s)
    du.x .= -1.0im * (hmat * u.x - u.x * hmat)
    ω_ba = QTBase.ω_matrix(D.H, D.lvl)
    D.Davies(du, u, ω_ba, p.tf, s)
end


function solve_af_rwa(
    A::Annealing,
    tf::Real;
    span_unit::Bool = false,
    ω_hint = [],
    lvl::Int = size(A.H, 1),
    kwargs...,
)
    if !(typeof(A.H) <: AdiabaticFrameHamiltonian)
        throw(ArgumentError("Adiabatic Frame RWA equation currently only works for adiabatic frame Hamiltonian."))
    end
    u0 = build_u0(A.u0, :m, control = A.control)
    tf = build_tf(tf, span_unit)
    #
    davies = build_davies(A.coupling, A.bath, ω_hint)
    f = AFRWADiffEqOperator(A.H, davies; lvl = lvl, control = A.control)
    p = LightAnnealingParams(tf; control = A.control)
    if typeof(A.control) <: PausingControl
        cb = DiscreteCallback(pause_condition, pause_affect!)
        kwargs = Dict{Symbol,Any}(kwargs)
        kwargs[:callback] = cb
    end
    tspan, tstops = scaling_time(tf, A.sspan, A.tstops)
    prob = ODEProblem{true}(f, u0, tspan, p)
    solve(prob; alg_hints = [:nonstiff], tstops = tstops, kwargs...)
end


function solve_af_rwa(
    A::Annealing,
    tf::Vector{T},
    alg,
    para_alg = EnsembleSerial();
    output_func = (sol, i) -> (sol, false),
    span_unit::Bool = false,
    ω_hint = [],
    lvl::Int = size(A.H, 1),
    kwargs...,
) where {T<:Real}
    if !(typeof(A.H) <: AdiabaticFrameHamiltonian)
        throw(ArgumentError("Adiabatic Frame RWA equation currently only works for adiabatic frame Hamiltonian."))
    end
    u0 = build_u0(A.u0, :m, control = A.control)
    t0 = build_tf(1.0, span_unit)
    davies = build_davies(A.coupling, A.bath, ω_hint)
    f = AFRWADiffEqOperator(A.H, davies; lvl = lvl, control = A.control)
    p = LightAnnealingParams(t0; control = A.control)
    # trajectories numbers
    trajectories = length(tf)
    tf_arr = float.(tf)
    # resolve control
    if typeof(A.control) <: PausingControl
        cb = DiscreteCallback(pause_condition, pause_affect!)
        kwargs = Dict{Symbol,Any}(kwargs)
        kwargs[:callback] = cb
    end
    #
    if span_unit == true
        tstops = create_tstops_for_tf_array(tf_arr, A.tstops)
        prob_func = (prob, i, repeat) -> begin
            tspan = (prob.tspan[1] * tf_arr[i], prob.tspan[2] * tf_arr[i])
            p = set_tf(prob.p, tf_arr[i])
            ODEProblem{true}(prob.f, prob.u0, tspan, p)
        end
    else
        tstops = A.tstops
        prob_func = (prob, i, repeat) -> begin
            p = set_tf(prob.p, tf_arr[i])
            ODEProblem{true}(prob.f, prob.u0, prob.tspan, p)
        end
    end
    prob = ODEProblem{true}(f, u0, A.sspan, p)
    ensemble_prob = EnsembleProblem(
        prob;
        prob_func = prob_func,
        output_func = output_func,
    )
    solve(
        ensemble_prob,
        alg,
        para_alg;
        trajectories = trajectories,
        tstops = tstops,
        kwargs...,
    )
end


function build_davies(coupling, bath::OhmicBath, ω_range, lambshift)
    γ_loc, S_loc = davies_spectrum(bath, ω_range, lambshift)
    DaviesGenerator(coupling, γ_loc, S_loc)
end


function build_ensemble_problem_ame_trajectory(
    A::Annealing,
    tf::Real,
    prob_func,
    output_func,
    reduction;
    span_unit::Bool = false,
    ω_hint = [],
    lambshift::Bool = true,
    lvl::Int = size(A.H, 1),
    eig_tol = 0.0,
    v0 = copy(A.u0),
    kwargs...,
)
    tf = build_tf(tf, span_unit)
    davies = build_davies(A.coupling, A.bath, ω_hint, lambshift)
    control = AMETrajectoryOperator(
        A.H,
        davies,
        lvl,
        eig_tol = eig_tol,
        v0 = v0,
    )
    control = A.control == nothing ? control : ControlSet(control, A.control)
    u0 = build_u0(A.u0, :v, control = control)
    p = LightAnnealingParams(tf; control = control)
    callback = build_callback(control, :ame_trajectory)
    ff = ame_trajectory_construct_ode_function(A.H, control)
    prob = ODEProblem{true}(ff, u0, (p) -> scaling_tspan(p.tf, A.sspan), p)

    ensemble_prob = EnsembleProblem(
        prob;
        prob_func = prob_func,
        output_func = output_func,
        reduction = reduction,
    )

    ensemble_prob, callback
end


function ame_trajectory_construct_ode_function(H, ::AMETrajectoryOperator)
    cache = get_cache(H)
    diff_op = DiffEqArrayOperator(
        cache,
        update_func = (A, u, p, t) -> update_cache!(A, p.control, p.tf, t),
    )
    jac_cache = similar(cache)
    jac_op = DiffEqArrayOperator(
        jac_cache,
        update_func = (A, u, p, t) -> update_cache!(A, p.control, p.tf, t),
    )
    ff = ODEFunction(diff_op; jac_prototype = jac_op)
end
