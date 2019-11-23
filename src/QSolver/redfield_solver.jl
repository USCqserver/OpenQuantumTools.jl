"""
    function solve_redfield(A::Annealing, tf::Real, unitary; span_unit::Bool = false, kwargs...)

Solve the time dependent Redfield equation for `Annealing` defined by `A` with total annealing time `tf`.

...
# Arguments
- `A::Annealing`: the Annealing object.
- `tf::Real`: the total annealing time.
- `unitary` : precalculated unitary of close system evolution.
- `span_unit::Bool=false`: flag variable which, when set to true, informs the solver to work with time in physical unit.
- `tstops` : extra times that the timestepping algorithm must step to.
- `kwargs` : other keyword arguments supported by DifferentialEquations.jl
...
"""
function solve_redfield(
    A::Annealing,
    tf::Real,
    unitary;
    span_unit::Bool = false,
    tstops = Float64[],
    positivity_check = false,
    kwargs...,
)
    tf = build_tf(tf, span_unit)
    tstops = build_tstops(tf, tstops, A.tstops)
    u0 = build_u0(A.u0, :m, control = A.control)
    ff = redfield_construct_ode_function(A.H, A.control)
    coupling = adjust_coupling_with_control(A.coupling, A.control)
    ff = redfield_construct_ode_function(A.H, A.control)
    opensys = create_redfield(coupling, unitary, tf, A.bath)
    p = AnnealingParams(A.H, tf; opensys = opensys, control = A.control)
    callback = build_callback(A.control, :redfield)
    if positivity_check
        positivity_check_callback = FunctionCallingCallback(
            positivity_check_affect,
            func_everystep = true,
            func_start = false,
        )
        callback = CallbackSet(callback, positivity_check_callback)
    end
    prob = ODEProblem(ff, u0, (p) -> scaling_tspan(p.tf, A.sspan), p)
    solve(prob; alg_hints = [:nonstiff], tstops = tstops, callback = callback, kwargs...)
end


function redfield_construct_ode_function(H, ::Union{Nothing,InstPulseControl})
    redfield_f
end


function redfield_contruct_ode_function(H, ::PausingControl)
    redfield_control_f
end


redfield_construct_coupling_function(coupling, ::Union{Nothing,InstPulseControl}) = coupling


redfield_construct_coupling_function(coupling, control::PausingControl) =
    attach_annealing_param(control, coupling)


function redfield_f(du, u, p, t)
    p.H(du, u, p.tf, t)
    p.opensys(du, u, p.tf, t)
end


function redfield_control_f(du, u, p, t)
    p.H(du, u, p, t)
    p.opensys(du, u, p.tf, t)
end
