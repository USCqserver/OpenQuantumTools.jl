@deprecate solve_davies(
    A::Annealing,
    tf::Real;
    span_unit = false, ω_hint = nothing, lvl = nothing, kwargs...
) solve_ame(A::Annealing, tf::Real; span_unit = false, ω_hint = nothing, lvl = nothing, kwargs...)

function solve_ame(
    A::Annealing,
    tf::Real;
    span_unit = false, ω_hint = nothing, lvl = nothing, kwargs...
)
    if ndims(A.u0) == 1
        u0 = A.u0 * A.u0'
    else
        u0 = A.u0
    end
    u0 = prepare_u0(u0, A.control)
    tf = prepare_tf(tf, span_unit)
    #
    davies = create_davies(A.coupling, A.bath; ω_range = ω_hint)
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


function create_davies(coupling, bath::OhmicBath; ω_range = nothing)
    γ_loc, S_loc = davies_spectrum(bath; ω_range = ω_range)
    DaviesGenerator(coupling, γ_loc, S_loc)
end


function (D::AMEDiffEqOperator{true,T})(du, u, p, t) where T <: PausingControl
    s, a_scale, g_scale = p.control(p.tf, t)
    hmat = D.H(u, a_scale, g_scale, s)
    ω_ba = ω_matrix(D.H)
    D.Davies(du, u, ω_ba, p.tf, t)
end
