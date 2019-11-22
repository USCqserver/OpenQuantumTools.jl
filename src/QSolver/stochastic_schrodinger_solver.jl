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
    tf = build_tf(tf, span_unit)
    tstops = build_tstops(tf, tstops, A.tstops)
    control = FluctuatorControl(tf, length(A.coupling), A.bath)
    control = merge_control(A.control, control)
    u0 = build_u0(A.u0, type = :v, control = control)
    opensys = create_redfield(coupling, unitary, tf, A.bath)
    p = AnnealingParams(A.H, tf; control = A.control)
    ff = schrodinger_construct_ode_function(A.H, A.control)
    prob = ODEProblem{true}(ff, u0, (p) -> scaling_tspan(p.tf, A.sspan), p)
    solve(prob; alg_hints = [:nonstiff], callback = callback, tstops = tstops, kwargs...)
end


function stochastic_schroding_construct_ode_function(H,)
end
