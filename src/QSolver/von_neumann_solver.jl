function solve_von_neumann(
    A::Annealing,
    tf::Real;
    dimensionless_time = true,
    vectorize::Bool = false,
    tstops = Float64[],
    de_array_constructor = nothing,
    kwargs...,
)
    tf, u0, tstops = __init(
        A,
        tf,
        dimensionless_time,
        :m,
        tstops,
        de_array_constructor,
        vectorize = vectorize,
    )
    reset!(A.control)
    callback = von_neumann_build_callback(A.control)
    ff = von_neumann_construct_ode_function(A.H, vectorize)
    p = ODEParams(A.H, tf; control = A.control)
    prob = ODEProblem{true}(ff, u0, (p) -> scaling_tspan(p.tf, A.sspan), p)
    solve(
        prob;
        alg_hints = [:nonstiff],
        tstops = tstops,
        callback = callback,
        kwargs...,
    )
end

von_neumann_build_callback(control) = nothing
von_neumann_build_callback(
    control::Union{InstPulseControl,InstDEPulseControl},
) = build_callback(control, pulse_on_density!)


function von_neumann_construct_ode_function(H, vectorize::Bool)
    cache = get_cache(H, vectorize)
    if vectorize == false
        ff = ODEFunction{true}(von_f; jac = von_jac)
    else
        diff_op = DiffEqArrayOperator(
            cache,
            update_func = (A, u, p, t) ->
                update_vectorized_cache!(A, p.H, p.tf, t),
        )
        jac_cache = similar(cache)
        jac_op = DiffEqArrayOperator(
            jac_cache,
            update_func = (A, u, p, t) ->
                update_vectorized_cache!(A, p.H, p.tf, t),
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
