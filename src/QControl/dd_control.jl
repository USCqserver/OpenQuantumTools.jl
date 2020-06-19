"""
$(TYPEDEF)
Controller object for applying instantaneous pulses during annealing.

# Fields
$(FIELDS)
"""
mutable struct InstPulseControl <: AbstractAnnealingControl
    """Positions of control pulses"""
    tstops::Any
    """Function to generate pulses based on its input"""
    pulse_func::Any
    """Index for internal state machine"""
    state::Int
end

InstPulseControl(tstops, pulse_func) = InstPulseControl(tstops, pulse_func, 1)
get_pulse(C::InstPulseControl) = C.pulse_func(C.state)
next!(C::InstPulseControl) = C.state += 1
reset!(C::InstPulseControl) = C.state = 1
reset!(C::InstPulseControl, u0, initializer) = C.state = 1
get_controller_name(::InstPulseControl) = :pulse_control


function build_callback(control::InstPulseControl, update!)
    affect! = function (integrator)
        controller = integrator.p.control
        pulse = get_pulse(controller)
        next!(controller)
        for c in full_cache(integrator)
            update!(c, pulse)
        end
    end
    PresetTimeCallback(control.tstops, affect!)
end


function build_callback(
    control::InstPulseControl,
    controller_name::Symbol,
    update!,
)
    affect! = function (integrator)
        controller = getproperty(integrator.p.control, controller_name)
        pulse = get_pulse(controller)
        next!(controller)
        for c in full_cache(integrator)
            update!(c, pulse)
        end
    end
    PresetTimeCallback(control.tstops, affect!)
end


"""
$(TYPEDEF)
Controller object for applying instantaneous pulses during annealing. This is to use with the `DEDataArray` structure to record the `state` information.

# Fields
$(FIELDS)
"""
struct InstDEPulseControl <: DEDataControl
    """Positions of control pulses"""
    tstops::Any
    """Function to generate pulses based on its input"""
    pulse_func::Any
    """Symbol used in DEDataArray"""
    sym::Symbol
end

get_pulse(control::InstDEPulseControl, state::Int) = control.pulse_func(state)
reset!(::InstDEPulseControl) = nothing
reset!(ctr::InstDEPulseControl, u0, initializer) = setfield!(u0, ctr.sym, 1)
get_controller_name(::InstDEPulseControl) = :pulse_control


function build_callback(control::InstDEPulseControl, update!)
    affect! = function (integrator)
        controller = integrator.p.control
        state_sym = controller.sym
        state = getproperty(integrator.u, state_sym)
        pulse = get_pulse(controller, state)
        for c in full_cache(integrator)
            update!(c.x, pulse)
            setproperty!(c, state_sym, state + 1)
        end
    end
    PresetTimeCallback(control.tstops, affect!)
end


function build_callback(
    control::InstDEPulseControl,
    controller_name::Symbol,
    update!,
)
    affect! = function (integrator)
        controller = getproperty(integrator.p.control, controller_name)
        state_sym = controller.sym
        state = getproperty(integrator.u, state_sym)
        pulse = get_pulse(controller, state)
        for c in full_cache(integrator)
            update!(c, pulse)
            setproperty!(c, state_sym, state + 1)
        end
    end
    PresetTimeCallback(control.tstops, affect!)
end
