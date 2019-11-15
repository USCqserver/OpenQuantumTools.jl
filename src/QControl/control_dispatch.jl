"""
    function construct_callback(control, solver_type)

Contruct the corresponding callback functions (if there is any) for `control`, which can be nothing or any `AbstractAnnealingControl` object. `solver_type` specifies for which solver the callback function is needed.
"""
function construct_callback(control::PausingControl, solver_type::Symbol)
    PresetTimeCallback(control.tstops, pause_affect!)
end


function construct_callback(control::InstPulseControl, solver_type::Symbol)
    if solver_type == :unitary
        PresetTimeCallback(control.tstops, unitary_dd_affect!)
    elseif solver_type == :redfield
        PresetTimeCallback(control.tstops, density_matrix_dd_affect!)
    else
        ArgumentError("$solver_type solver does not support the specified control protocol.")
    end
end


function construct_callback(control::FluctuatorControl, solver_type::Symbol)
    if solver_type == :stochastic_schrodinger
        IterativeCallback()
    else
        ArgumentError("$solver_type solver does not support the specified control protocol.")
    end
end

"""
    function adjust_u0_with_control(u0, p)

Convert the state vector/density matrix to the corresponding DEDataArray depending on the type of control `p`.
"""
function adjust_u0_with_control(u0, ::Union{PausingControl,InstPulseControl})
    if ndims(u0) == 1
        DEStateMachineVec(u0, 1)
    elseif ndims(u0) == 2
        DEStateMachineMat(u0, 1)
    else
        throw(ArgumentError("u0 can either be a vector or matrix."))
    end
end


"""
    function pause_affect!(integrator)

Callback function for `PausingControl`.
"""
function pause_affect!(integrator)
    for c in full_cache(integrator)
        c.state = 1 - c.state
    end
    u_modified!(integrator, false)
end


"""
    function unitary_dd_affect!(integrator)

Callback function for `InstPulseControl`. It is used by unitary solver.
"""
function unitary_dd_affect!(integrator)
    for c in full_cache(integrator)
        pulse = integrator.p.control(c.state)
        c.x = pulse_on_unitary(pulse, c)
        c.state += 1
    end
end


@inline pulse_on_unitary(p, c::DEDataVector) = Matrix{eltype(p)}(I, size(p))⊗p * c.x
@inline pulse_on_unitary(p, c::DEDataMatrix) = p * c.x


function density_matrix_dd_affect!(integrator)
    for c in full_cache(integrator)
        pulse = integrator.p.control(c.state)
        c.x = pulse_on_density_matrix(pulse, c)
        c.state += 1
    end
end

@inline pulse_on_density_matrix(p, c::DEDataVector) = conj(p)⊗p * c.x
@inline pulse_on_density_matrix(p, c::DEDataMatrix) = p * c.x * p'


function fluctuator_affect!(integrator)
    noise_value = integrator.p.control()
    for c in full_cache(integrator)
        c.n = noise_value
    end
    next_state!(integrator.p.control)
    u_modified!(integrator, false)
end


function fluctuator_time_choice(integrator)
    integrator.t + integrator.p.control.next_τ
end
