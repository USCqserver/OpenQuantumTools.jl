function solve_stochastic_schrodinger(
    A::Annealing,
    tf::Real,
    alg,
    para_alg = EnsembleSerial();
    output_func = (sol, i) -> (sol, false),
    span_unit = false,
    tstops = Float64[],
    kwargs...,
)
    tf = prepare_tf(tf, span_unit)
    tstops = prepare_tstops(tf, tstops, A.tstops)
    control = construct_stochastic_control(tf, A.bath)
    opensys = create_redfield(coupling, unitary, tf, A.bath)
    p = AnnealingParams(A.H, tf; control = A.control)
    callback = construct_callback(A.control, :schrodinger)
    ff = schrodinger_construct_ode_function(A.H, A.control)
    prob = ODEProblem{true}(ff, u0, (p) -> scaling_tspan(p.tf, A.sspan), p)
    solve(prob; alg_hints = [:nonstiff], callback = callback, tstops = tstops, kwargs...)
end


function stochastic_schroding_construct_ode_function(H,)
end
