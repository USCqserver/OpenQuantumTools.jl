function solve_unitary(A::Annealing, tf::Real; span_unit = false, kwargs...)
    tstops = A.tstops
    u0 = Matrix{ComplexF64}(I, A.H.size)
    u0 = prepare_u0(u0, A.control)
    tf = prepare_tf(tf, span_unit)
    #jp = vectorized_jacobian_prototype(A.H)
    p = AnnealingParams(A.H, tf; control = A.control)
    if typeof(A.control) <: PausingControl
        ff = ODEFunction(uni_control_f; jac = uni_control_jac)
        cb = DiscreteCallback(pause_condition, pause_affect!)
        kwargs = Dict{Symbol,Any}(kwargs)
        kwargs[:callback] = cb
    else
        ff = ODEFunction(uni_f; jac = uni_jac)
    end
    tspan, tstops = scaling_time(tf, A.sspan, tstops)
    prob = ODEProblem{true}(ff, u0, tspan, p)
    solve(prob; alg_hints = [:nonstiff], tstops = tstops, kwargs...)
end


function uni_f(du, u, p, t)
    hmat = p.H(p.tf, t)
    mul!(du, hmat, u)
    lmul!(-1.0im, du)
end


function uni_jac(J, u, p, t)
    hmat = p.H(p.tf, t)
    mul!(J, -1.0im, Matrix{eltype(hmat)}(I, p.H.size) ⊗ hmat)
end


function uni_control_f(du, u, p::AbstractAnnealingParams, t::Real)
    hmat = p.H(u, p, t)
    mul!(du, hmat, u)
    lmul!(-1.0im, du)
end


function uni_control_jac(J, u, p::AbstractAnnealingParams, t::Real)
    hmat = p.H(u, p, t)
    mul!(J, -1.0im, Matrix{eltype(hmat)}(I, p.H.size) ⊗ hmat)
end
