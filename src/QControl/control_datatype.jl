abstract type DEDataControl <: AbstractAnnealingControl end
get_symbol(control::DEDataControl) = control.sym

require_de_data(::DEDataControl) = true
require_de_data(::AbstractAnnealingControl) = false
require_de_data(::Nothing) = false

(h::DenseHamiltonian)(du::DEDataArray, u::DEDataArray, tf, t) =
    (h::DenseHamiltonian)(du.x, u.x, tf, t)
QTBase.check_positivity(data::DEDataArray) = check_positivity(data.x)
reset!(::Nothing) = nothing
struct DEFAULT_INITIALIZER_CLASS end
const DEFAULT_INITIALIZER = DEFAULT_INITIALIZER_CLASS()


"""
$(TYPEDEF)
ControlSet to combine multiple control protocols.

# Fields
$(FIELDS)
"""
struct ControlSet{T<:NamedTuple} <: AbstractAnnealingControl
    """NamedTuple for controllers"""
    ctrs::T
end


function ControlSet(names, controllers)
    ControlSet((; zip(names, controllers)...))
end


function ControlSet(args::AbstractAnnealingControl...)
    names = [get_controller_name(c) for c in args]
    if !allunique(names)
        error("Using the multiple controllers of the same type is currently unsupported.")
    end
    ControlSet(names, args)
end


function ControlSet(c::ControlSet, args::AbstractAnnealingControl...)
    ctrs = getfield(c, :ctrs)
    set_names = [v for v in keys(ctrs)]
    names = vcat(set_names, [get_controller_name(c) for c in args])

    if !allunique(names)
        error("Using the multiple controllers of the same type is currently unsupported.")
    end
    controllers = vcat(ctrs..., args...)

    ControlSet(names, controllers)
end


ControlSet(::Nothing) = nothing

Base.getproperty(x::ControlSet, sym::Symbol) = getfield(getfield(x, :ctrs), sym)
Base.length(c::ControlSet) = length(getfield(c, :ctrs))
Base.iterate(c::ControlSet, state = 1) = Base.iterate(getfield(c, :ctrs), state)
Base.keys(c::ControlSet) = Base.keys(getfield(c, :ctrs))

reset!(ctrs::ControlSet, u0, initializer) =
    [reset!(c, u0, initializer) for c in ctrs]

has_fluctuator_control(c::ControlSet) =
    haskey(getfield(c, :ctrs), :fluctuator_control)
has_ame_trajectory_control(c::ControlSet) =
    haskey(getfield(c, :ctrs), :ame_trajectory_control)
need_de_array(c::ControlSet, sym::Symbol) =
    typeof(getproperty(c, sym)) <: DEDataControl