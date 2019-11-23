"""
    function build_callback(control, solver_type)

Contruct the corresponding callback functions (if there is any) for `control`, which can be nothing or any `AbstractAnnealingControl` object. `solver_type` specifies for which solver the callback function is needed.
"""
function build_callback(::Nothing, ::Symbol)
    nothing
end


function build_callback(control::PausingControl, solver_type::Symbol)
    PresetTimeCallback(control.tstops, pause_affect!)
end


function build_callback(control::InstPulseControl, solver_type::Symbol)
    if solver_type == :unitary
        PresetTimeCallback(control.tstops, unitary_dd_affect!)
    elseif solver_type == :redfield
        PresetTimeCallback(control.tstops, density_matrix_dd_affect!)
    else
        ArgumentError("$solver_type solver does not support the specified control protocol.")
    end
end


function build_callback(control::FluctuatorControl, solver_type::Symbol)
    if solver_type == :stochastic_schrodinger
        IterativeCallback(fluctuator_time_choice, fluctuator_affect!)
    else
        ArgumentError("$solver_type solver does not support the specified control protocol.")
    end
end


function build_callback(control::ControlSet, solver_type::Symbol)
    if length(control) == 1
        res = build_callback(control.ctrs[1], solver_type)
    else
        cbs = [build_callback(c, solver_type) for c in control]
        res = CallbackSet(cbs...)
    end
    res
end


"""
    function adjust_u0_with_control(u0, p)

Convert the state vector/density matrix to the corresponding DEDataArray depending on the type of control `p`.
"""
function adjust_u0_with_control(u0, ::Union{Nothing, ParameterFreeControl})
    u0
end


function adjust_u0_with_control(
    u0::Array{T,N},
    ::Union{PausingControl,InstPulseControl},
) where {T<:Number,N}
    DEStateMachineArray{T,N}(u0, 1)
end


function adjust_u0_with_control(u0::Array{T,N}, f::FluctuatorControl) where {T<:Number,N}
    DENoiseArray{T,N}(u0, f())
end


function adjust_u0_with_control(u0, control::ControlSet)
    if length(control)==1
        adjust_u0_with_control(u0, control.ctrs[1])
    else
        error("Adjusting initial state according to a ControlSet is currently not implemented.")
    end
end


function adjust_coupling_with_control(coupling, ::Union{Nothing, InstPulseControl})
    coupling
end


function adjust_coupling_with_control(coupling, control::PausingControl)
    attach_annealing_param(control, coupling)
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
