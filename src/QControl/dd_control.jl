"""
$(TYPEDEF)
Controller object for applying instantaneous pulses during annealing

# Fields
$(FIELDS)
"""
struct InstPulseControl <: AbstractAnnealingControl
    """Positions of control pulses"""
    tstops
    """Function to generate pulses based on its input"""
    pulse_func
end


function adjust_u0_with_control(u0, ::InstPulseControl)
    if ndims(u0) == 1
        DEStateMachineVec(u0, 1)
    elseif ndims(u0) == 2
        DEStateMachineMat(u0, 1)
    else
        throw(ArgumentError("u0 can either be a vector or matrix."))
    end
end



function create_callback(control::InstPulseControl, update_func)
    user_affect! =
    PresetTimeCallback(pulse_position, )
end


function pulse_unitary_affect!(integrator, p_list)
  for c in full_cache(integrator)
    c.state += 1
    lmul!(p_list[c.state], c.x)
  end
end
