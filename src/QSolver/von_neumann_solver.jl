function solve_von_neumann(
    A::Annealing,
    tf::Real;
    span_unit = false,
    vectorize = false,
    tstops = Float64[],
    kwargs...,
)
    tf = build_tf(tf, span_unit)
    tstops = build_tstops(tf, tstops, A.tstops)
    u0 = build_u0(A.u0, :m, control = A.control, vectorize = vectorize)
    callback = build_callback(A.control, :von_neumann)
    ff = von_neumann_construct_ode_function(A.H, A.control, vectorize)
    p = AnnealingParams(A.H, tf; control = A.control)
    prob = ODEProblem{true}(ff, u0, (p) -> scaling_tspan(p.tf, A.sspan), p)
    solve(prob; alg_hints = [:nonstiff], tstops = tstops, callback = callback, kwargs...)
end


function von_neumann_construct_ode_function(
    H,
    control::Union{Nothing,InstPulseControl},
    vectorize::Bool,
)
    cache = get_cache(H, vectorize)
    if vectorize == false
        ff = ODEFunction{true}(von_f; jac = von_jac)
    else
        diff_op = DiffEqArrayOperator(
            cache,
            update_func = (A, u, p, t) -> update_vectorized_cache!(A, p.H, p.tf, t),
        )
        jac_cache = similar(cache)
        jac_op = DiffEqArrayOperator(
            jac_cache,
            update_func = (A, u, p, t) -> update_vectorized_cache!(A, p.H, p.tf, t),
        )
        ff = ODEFunction(diff_op; jac_prototype = jac_op)
    end
    ff
end


function von_f(du, u, p, t)
    p.H(du, u, p.tf, t)
end


function von_jac(J, u, p, t)
    hmat = p.H(p.tf, t)
    iden = Matrix{eltype(hmat)}(I, p.H.size)
    temp = (iden ⊗ hmat - transpose(hmat) ⊗ iden)
    mul!(J, -1.0im, temp)
end

function von_control_f(du, u, p, t)
    p.H(du, u, p, t)
end

function von_control_jac(J, u, p, t)
    hmat = p.H(u, p, t)
    iden = Matrix{eltype(hmat)}(I, p.H.size)
    temp = (iden ⊗ hmat - transpose(hmat) ⊗ iden)
    mul!(J, -1.0im, temp)
end
