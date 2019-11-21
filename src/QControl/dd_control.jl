"""
$(TYPEDEF)
Controller object for applying instantaneous pulses during annealing.

# Fields
$(FIELDS)
"""
struct InstPulseControl <: AbstractAnnealingControl
    """Positions of control pulses"""
    tstops
    """Function to generate pulses based on its input"""
    pulse_func
end

(C::InstPulseControl)(state::Int) = C.pulse_func(state)

QTBase.need_change_time_scale(::InstPulseControl) = false
