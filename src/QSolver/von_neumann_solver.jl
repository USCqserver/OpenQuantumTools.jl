function solve_von_neumann(A::Annealing, tf::Real; span_unit = false, kwargs...)
    if ndims(A.u0) == 1
        u0 = A.u0*A.u0'
    else
        u0 = A.u0
    end
    tstops = A.tstops
    u0 = prepare_u0(u0, span_unit)
    tf = prepare_tf(tf, span_unit)
    jp = prepare_jacobian_prototype(A.H)
    p = AnnealingParams(A.H, tf; control = A.control)
    if typeof(A.control) <: PausingControl
        ff = ODEFunction(
            von_control_f;
            jac = von_control_jac, jac_prototype = jp
        )
        cb = DiscreteCallback(pause_condition, pause_affect!)
        kwargs = Dict{Symbol, Any}(kwargs)
        kwargs[:callback] = cb
    else
        ff = ODEFunction(von_f; jac = von_jac, jac_prototype = jp)
    end
    tspan, tstops = scaling_time(tf, A.sspan, tstops)
    prob = ODEProblem{true}(ff, u0, tspan, p)
    solve(prob; alg_hints = [:nonstiff], tstops = tstops, kwargs...)
end

function von_f(du, u, p, t)
end

function von_control_f()
end
