"""
    function solve_ame(
        A::Annealing,
        tf::Real;
        span_unit::Bool = false,
        ω_hint = [],
        lvl = size(A.H, 1),
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
- `lvl::Int = size(A.H, 1)` : number of levels to keep.
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
    kwargs...,
)
    u0 = build_u0(A.u0, :m, control = A.control)
    tf = build_tf(tf, span_unit)
    davies = build_davies(A.coupling, A.bath, ω_hint, lambshift)
    f = AMEDiffEqOperator(A.H, davies; lvl = lvl, control = A.control)
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


function (D::AMEDiffEqOperator{true,T})(du, u, p, t) where {T<:PausingControl}
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
    ensemble_prob = EnsembleProblem(prob; prob_func = prob_func, output_func = output_func)
    solve(ensemble_prob, alg, para_alg; trajectories = trajectories, tstops = tstops, kwargs...)
end


function build_davies(coupling, bath::OhmicBath, ω_range, lambshift)
    γ_loc, S_loc = davies_spectrum(bath, ω_range, lambshift)
    DaviesGenerator(coupling, γ_loc, S_loc)
end
