"""
    function solve_redfield(
        A::Annealing,
        tf::Real,
        unitary;
        vectorize::Bool = false,
        dimensionless_time::Bool = true,
        tstops = Float64[],
        positivity_check::Bool = false,
        de_array_constructor = nothing,
        int_atol = 1e-8,
        int_rtol = 1e-6,
        kwargs...,
    )

Solve the time dependent Redfield equation for `Annealing` defined by `A` with total annealing time `tf`.

...
# Arguments
- `A::Annealing`: the Annealing object.
- `tf::Real`: the total annealing time.
- `unitary`: precalculated unitary of close system evolution.
- `vectorize::Bool = false`: whether to vectorize the density matrix.
- `dimensionless_time::Bool=true`: flag variable which, when set to true, informs the solver to work with dimensionless time.
- `tstops`: extra times that the timestepping algorithm must step to.
- `positivity_check::Bool = false`: whether to check the positivity of density matrix at each time step.
- `de_array_constructor = nothing`: the converting function if using `DEDataArray` type.
- `int_atol = 1e-8`: the absolute error tolerance for integration.
- `int_rtol = 1e-6`: the relative error tolerance for integration.
- `kwargs`: other keyword arguments supported by DifferentialEquations.jl.
...
"""
function solve_redfield(
    A::Annealing,
    tf::Real,
    unitary;
    vectorize::Bool = false,
    dimensionless_time::Bool = true,
    tstops = Float64[],
    positivity_check::Bool = false,
    de_array_constructor = nothing,
    int_atol = 1e-8,
    int_rtol = 1e-6,
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
    #coupling = adjust_coupling_with_control(A.coupling, A.control)
    ff = redfield_construct_ode_function(A.H, vectorize)
    opensys = build_redfield(
        A.interactions,
        unitary,
        tf,
        atol = int_atol,
        rtol = int_rtol,
    )
    reset!(A.control)
    callback = redfield_build_callback(A.control)
    p = ODEParams(A.H, tf; opensys = opensys, control = A.control)
    if positivity_check
        positivity_check_callback = FunctionCallingCallback(
            positivity_check_affect,
            func_everystep = true,
            func_start = false,
        )
        callback = CallbackSet(callback, positivity_check_callback)
    end
    prob = ODEProblem(ff, u0, (p) -> scaling_tspan(p.tf, A.sspan), p)
    solve(
        prob;
        alg_hints = [:nonstiff],
        tstops = tstops,
        callback = callback,
        kwargs...,
    )
end


redfield_build_callback(control) = nothing
redfield_build_callback(control::Union{InstPulseControl,InstDEPulseControl}) =
    build_callback(control, pulse_on_density!)

function redfield_construct_ode_function(H, vectorize)
    if vectorize == false
        redfield_f
    else
        cache = get_cache(H, vectorize)
        diff_op =
            DiffEqArrayOperator(cache, update_func = redfield_vectorize_update!)
        jac_cache = similar(cache)
        jac_op = DiffEqArrayOperator(
            jac_cache,
            update_func = redfield_vectorize_update!,
        )
        ff = ODEFunction(diff_op; jac_prototype = jac_op)
    end
end


function redfield_f(du, u, p, t)
    p.H(du, u, p.tf, t)
    p.opensys(du, u, p.tf, t)
end


function redfield_vectorize_update!(A, u, p, t)
    update_vectorized_cache!(A, p.H, p.tf, t)
    update_vectorized_cache!(A, p.opensys, p.tf, t)
end
