function solve_schrodinger(
    A::Annealing,
    tf::Real;
    dimensionless_time = true,
    manifold_projection = false,
    tstops = Float64[],
    de_array_constructor = nothing,
    kwargs...,
)
    tf, u0, tstops = __init(
        A,
        tf,
        dimensionless_time,
        :v,
        tstops,
        de_array_constructor,
    )
    reset!(A.control)
    callback = schrodinger_build_callback(A.control)
    p = ODEParams(A.H, tf; control = A.control)

    if manifold_projection
        cb = ManifoldRetraction(Sphere())
        callback = callback == nothing ? cb : CallbackSet(callback, cb)
    end

    cache = get_cache(A.H)
    diff_op = DiffEqArrayOperator(
        cache,
        update_func = (A, u, p, t) -> update_cache!(A, p.H, p.tf, t),
    )
    jac_cache = similar(cache)
    jac_op = DiffEqArrayOperator(
        jac_cache,
        update_func = (A, u, p, t) -> update_cache!(A, p.H, p.tf, t),
    )
    ff = ODEFunction(diff_op; jac_prototype = jac_op)

    prob = ODEProblem{true}(ff, u0, (p) -> scaling_tspan(p.tf, A.sspan), p)
    solve(
        prob;
        alg_hints = [:nonstiff],
        callback = callback,
        tstops = tstops,
        kwargs...,
    )
end

schrodinger_build_callback(control) = nothing
schrodinger_build_callback(
    control::Union{InstPulseControl,InstDEPulseControl},
) = build_callback(control, (c, pulse) -> c .= pulse * c)
schrodinger_build_callback(control::ControlSet) =
    error("Control set is not currently supported for Schrodinger solver.")


# function schrodinger_construct_ode_function(
#     H,
#     ::Union{PausingControl,PausingDEControl},
# )
#     cache = get_cache(H)
#     diff_op = DiffEqArrayOperator(
#         cache,
#         update_func = (A, u, p, t) -> update_cache!(A, p.H, p.tf, t),
#     )
#     jac_cache = similar(cache)
#     jac_op = DiffEqArrayOperator(
#         jac_cache,
#         update_func = (A, u, p, t) -> update_cache!(A, p.H, p.tf, t),
#     )
#     ff = ODEFunction(diff_op; jac_prototype = jac_op)
# end
#
#
# function schrodinger_construct_ode_function(H, pause_control::PausingControl)
#     cache = get_cache(H)
#     diff_op = DiffEqArrayOperator(
#         cache,
#         update_func = (A, u, p, t) ->
#             update_cache!(A, u, p.tf, t, p.H, p.control),
#     )
#     jac_cache = similar(cache)
#     jac_op = DiffEqArrayOperator(
#         jac_cache,
#         update_func = (A, u, p, t) ->
#             update_cache!(A, u, p.tf, t, p.H, p.control),
#     )
#     ff = ODEFunction(diff_op; jac_prototype = jac_op)
# end
