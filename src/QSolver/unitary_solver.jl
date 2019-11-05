function solve_unitary(
    A::Annealing,
    tf::Real;
    span_unit = false,
    vectorize = false,
    kwargs...,
)
    tstops = A.tstops
    u0 = Matrix{ComplexF64}(I, A.H.size)
    u0 = prepare_u0(u0, type = :m, control = A.control, vectorize=vectorize)
    tf = prepare_tf(tf, span_unit)
    #jp = vectorized_jacobian_prototype(A.H)
    p = AnnealingParams(A.H, tf; control = A.control)
    if typeof(A.control) <: PausingControl
        ff = ODEFunction(uni_control_f; jac = uni_control_jac)
        cb = DiscreteCallback(pause_condition, pause_affect!)
        kwargs = Dict{Symbol,Any}(kwargs)
        kwargs[:callback] = cb
    else
        ff = uni_create_ode_fun(A.H, vectorize)
    end
    tspan, tstops = scaling_time(tf, A.sspan, tstops)
    prob = ODEProblem{true}(ff, u0, tspan, p)
    solve(prob; alg_hints = [:nonstiff], tstops = tstops, kwargs...)
end


function uni_create_ode_fun(H, vectorize)
    j_cache = Matrix{eltype(H)}(I, size(H.u_cache)) ⊗ H.u_cache
    if vectorize == false
        cache = H.u_cache
        diff_op_update = (A, u, p, t) -> update_cache!(A, p.H, p.tf, t)
    else
        cache = Matrix{eltype(H)}(I, size(H.u_cache)) ⊗ H.u_cache
        diff_op_update = uni_jac!
    end
    diff_op = DiffEqArrayOperator(cache, update_func = diff_op_update)
    jac_op = DiffEqArrayOperator(
        j_cache,
        update_func = uni_jac!,
    )
    ff = ODEFunction(diff_op; jac_prototype = jac_op)
end


function uni_jac!(J, u, p, t)
    hmat = p.H(p.tf, t)
    J .= -1.0im * Matrix{eltype(hmat)}(I, p.H.size) ⊗ hmat
end


#TODO Replace the code below with DiffEqOperators
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
