"""
    function solve_redfield(A::Annealing, tf::Real, unitary; span_unit = false, kwargs...)

Solve the time dependent Redfield equation for `Annealing` defined by `A` with total annealing time `tf`. The third argument is the precalculated unitary of close system evolution. It needs to be a callable object which returns the value of unitary at given time. `span_unit` is a boolean flag which, when set to true, informs the solver to work with time in physical unit.
"""
function solve_redfield(A::Annealing, tf::Real, unitary; span_unit = false, kwargs...)
    u0 = prepare_u0(A.u0, type=:m, control=A.control)
    tf = prepare_tf(tf, span_unit)
    coupling = A.coupling
    ff = redfield_f
    if typeof(A.control) <: PausingControl
        cb = DiscreteCallback(pause_condition, pause_affect!)
        kwargs = Dict{Symbol,Any}(kwargs)
        kwargs[:callback] = cb
        coupling = attach_annealing_param(A.control, coupling)
        ff = redfield_control_f
    end
    opensys = create_redfield(coupling, unitary, tf, A.bath)
    p = AnnealingParams(A.H, tf; opensys=opensys, control=A.control)
    tspan, tstops = scaling_time(tf, A.sspan, A.tstops)
    prob = ODEProblem(ff, u0, tspan, p)
    solve(prob; alg_hints = [:nonstiff], tstops=tstops, kwargs...)
end


function redfield_f(du, u, p, t)
    p.H(du, u, p.tf, t)
    p.opensys(du, u, p.tf, t)
end


function redfield_control_f(du, u, p, t)
    p.H(du, u, p, t)
    p.opensys(du, u, p.tf, t)
end
