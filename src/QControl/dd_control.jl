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
reset!(C::InstPulseControl) = nothing
QTBase.need_change_time_scale(::InstPulseControl) = false


"""
    function unitary_dd_affect!(integrator)

Callback function for `InstPulseControl`. It is used by unitary solver.
"""
function unitary_dd_affect!(integrator)
    pulse = integrator.p.control(integrator.u.state)
    for c in full_cache(integrator)
        c.x = pulse_on_unitary(pulse, c)
        c.state += 1
    end
end


@inline pulse_on_unitary(p, c::DEDataVector) = Matrix{eltype(p)}(I, size(p)) ⊗ p * c.x
@inline pulse_on_unitary(p, c::DEDataMatrix) = p * c.x


function density_matrix_dd_affect!(integrator)
    pulse = integrator.p.control(integrator.u.state)
    for c in full_cache(integrator)
        c.x = pulse_on_density_matrix(pulse, c)
        c.state += 1
    end
end

@inline pulse_on_density_matrix(p, c::DEDataVector) = conj(p) ⊗ p * c.x
@inline pulse_on_density_matrix(p, c::DEDataMatrix) = p * c.x * p'
