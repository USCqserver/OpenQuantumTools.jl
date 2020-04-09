"""
    function solve_ame(
        A::Annealing,
        tf::Real;
        dimensionless_time::Bool = true,
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
- `dimensionless_time::Bool=true`: flag variable which, when set to true, informs the solver to work with dimensionless time.
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
    dimensionless_time::Bool = true,
    ω_hint = [],
    lambshift::Bool = true,
    lvl::Int = size(A.H, 1),
    eig_tol = 0.0,
    maxiter::Int = 3000,
    ncv = 20,
    v0 = copy(A.u0),
    kwargs...,
)
    tf, tstops = preprocessing_time(tf, tstops, A.tstops, dimensionless_time)
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
    p = ODEParams(tf; control = A.control)
    # if typeof(A.control) <: PausingControl
    #     cb = DiscreteCallback(pause_condition, pause_affect!)
    #     kwargs = Dict{Symbol,Any}(kwargs)
    #     kwargs[:callback] = cb
    # end
    prob = ODEProblem(f, u0, (p) -> scaling_tspan(p.tf, A.sspan), p)
    solve(prob; alg_hints = [:nonstiff], tstops = tstops, kwargs...)
end


# function (D::QTBase.AMEDenseDiffEqOperator{true,T})(
#     du,
#     u,
#     p,
#     t,
# ) where {T<:PausingControl}
#     s, a_scale, g_scale = p.control(p.tf, t)
#     hmat = D.H(u, a_scale, g_scale, s)
#     du.x .= -1.0im * (hmat * u.x - u.x * hmat)
#     ω_ba = QTBase.ω_matrix(D.H, D.lvl)
#     D.Davies(du, u, ω_ba, p.tf, s)
# end


# function solve_af_rwa(
#     A::Annealing,
#     tf::Real;
#     dimensionless_time::Bool = true,
#     ω_hint = [],
#     lvl::Int = size(A.H, 1),
#     kwargs...,
# )
#     if !(typeof(A.H) <: AdiabaticFrameHamiltonian)
#         throw(ArgumentError("Adiabatic Frame RWA equation currently only works for adiabatic frame Hamiltonian."))
#     end
#     tf, tstops = preprocessing_time(tf, tstops, A.tstops, dimensionless_time)
#     u0 = build_u0(A.u0, :m, control = A.control)
#     #
#     davies = build_davies(A.coupling, A.bath, ω_hint)
#     f = AFRWADiffEqOperator(A.H, davies; lvl = lvl, control = A.control)
#     p = LightAnnealingParams(tf; control = A.control)
#     if typeof(A.control) <: PausingControl
#         cb = DiscreteCallback(pause_condition, pause_affect!)
#         kwargs = Dict{Symbol,Any}(kwargs)
#         kwargs[:callback] = cb
#     end
#     prob = ODEProblem{true}(f, u0, (p) -> scaling_tspan(p.tf, A.sspan), p)
#     solve(prob; alg_hints = [:nonstiff], tstops = tstops, kwargs...)
# end


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
    dimensionless_time::Bool = true,
    ω_hint = [],
    lambshift::Bool = true,
    lvl::Int = size(A.H, 1),
    eig_tol = 0.0,
    maxiter::Int = 3000,
    ncv = 20,
    v0 = copy(A.u0),
    kwargs...,
)
    tf, tstops = preprocessing_time(tf, tstops, A.tstops, dimensionless_time)
    davies = build_davies(A.coupling, A.bath, ω_hint, lambshift)
    control =
        AMETrajectoryOperator(A.H, davies, lvl, eig_tol = eig_tol, v0 = v0)
    control = A.control == nothing ? control : ControlSet(control, A.control)
    u0 = build_u0(A.u0, :v, control = control)
    p = ODEParams(tf; control = control)
    callback = ame_trajectory_build_callback(control)
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


function ame_trajectory_build_callback(control)
    AMEJumpCallback()
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
