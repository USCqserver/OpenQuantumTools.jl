function solve_unitary(
    A::Annealing,
    tf::Real;
    span_unit = false,
    vectorize = false,
    tstops = Float64[],
    kwargs...,
)
    tf = prepare_tf(tf, span_unit)
    tstops = prepare_tstops(tf, tstops, A.tstops)
    u0 = Matrix{ComplexF64}(I, size(A.H))
    u0 = prepare_u0(u0, type = :m, control = A.control, vectorize = vectorize)
    p = AnnealingParams(A.H, tf; control = A.control)
    callback = construct_callback(A.control, :unitary)
    ff = unitary_construct_ode_function(A.H, A.control, vectorize)
    prob = ODEProblem{true}(ff, u0, (p) -> scaling_tspan(p.tf, A.sspan), p)
    solve(prob; alg_hints = [:nonstiff], tstops = tstops, callback = callback, kwargs...)
end


"""
    function unitary_construct_ode_function(H, control, vectorize)

Construct the corresponding `ODEFunction` for unitary solver.
"""
function unitary_construct_ode_function(H, ::Union{Nothing,InstPulseControl}, vectorize)
    cache = get_cache(H)
    j_cache = Matrix{eltype(H)}(I, size(H)) ⊗ cache
    if vectorize == false
        diff_op_update = (A, u, p, t) -> update_cache!(A, p.H, p.tf, t)
    else
        cache = Matrix{eltype(H)}(I, size(H)) ⊗ cache
        diff_op_update = uni_jac!
    end
    diff_op = DiffEqArrayOperator(cache, update_func = diff_op_update)
    jac_op = DiffEqArrayOperator(j_cache, update_func = uni_jac!)
    ff = ODEFunction(diff_op; jac_prototype = jac_op)
end


function unitary_construct_ode_function(H, control::PausingControl, vectorize)
    if !(typeof(H) <: AdiabaticFrameHamiltonian)
        throw(ArgumentError("Pausing control is only for adiabatic frame Hamiltonian. Use `tstops` for other type of Hamiltonians."))
    end
    cache = get_cache(H)
    j_cache = Matrix{eltype(H)}(I, size(H)) ⊗ cache
    if vectorize == false
        diff_op_update = (A, u, p, t) -> update_cache!(A, u, p.tf, t, p.H, p.control)
    else
        cache = Matrix{eltype(H)}(I, size(H)) ⊗ cache
        diff_op_update = uni_control_jac!
    end
    diff_op = DiffEqArrayOperator(cache, update_func = diff_op_update)
    jac_op = DiffEqArrayOperator(j_cache, update_func = uni_control_jac!)
    ff = ODEFunction(diff_op; jac_prototype = jac_op)
end


"""
    function uni_jac!(J, u, p, t)

Update the jacobian matrix for unitary solver. The Jacobian matrix should be able to multiply the vectorized version of the unitary.
"""
function uni_jac!(J, u, p, t)
    hmat = p.H(p.tf, t)
    J .= -1.0im * Matrix{eltype(hmat)}(I, size(p.H)) ⊗ hmat
end


function uni_control_jac!(J, u, p::AbstractAnnealingParams, t::Real)
    hmat = p.H(u, tf, t, p.control)
    J .= -1.0im * Matrix{eltype(hmat)}(I, size(p.H)) ⊗ hmat
end
