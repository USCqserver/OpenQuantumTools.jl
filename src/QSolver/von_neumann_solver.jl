function solve_von_neumann(A::Annealing, tf::Real; span_unit = false, kwargs...)
    tstops = A.tstops
    if ndims(A.u0) == 1
        u0 = A.u0*A.u0'
    else
        u0 = A.u0
    end
    u0 = prepare_u0(u0, A.control)
    tf = prepare_tf(tf, span_unit)
    if typeof(A.control) <: PausingControl
        ff = ODEFunction{true}(
            von_control_f;
            jac = von_control_jac
        )
        cb = DiscreteCallback(pause_condition, pause_affect!)
        kwargs = Dict{Symbol, Any}(kwargs)
        kwargs[:callback] = cb
    else
        ff = ODEFunction{true}(von_f; jac = von_jac)
    end
    tspan, tstops = scaling_time(tf, A.sspan, tstops)
    p = AnnealingParams(A.H, tf; control = A.control)
    prob = ODEProblem{true}(ff, u0, tspan, p)
    solve(prob; alg_hints = [:nonstiff], tstops = tstops, kwargs...)
end


function von_f(du, u, p, t)
    p.H(du, u, p.tf, t)
end


function von_jac(J, u, p, t)
    hmat = p.H(p.tf, t)
    iden = Matrix{eltype(hmat)}(I, p.H.size)
    temp = (iden⊗hmat - transpose(hmat)⊗iden)
    mul!(J, -1.0im, temp)
end

function von_control_f(du, u, p, t)
    p.H(du, u, p, t)
end

function von_control_jac(J, u, p, t)
    hmat = p.H(u, p, t)
    iden = Matrix{eltype(hmat)}(I, p.H.size)
    temp = (iden⊗hmat - transpose(hmat)⊗iden)
    mul!(J, -1.0im, temp)
end
