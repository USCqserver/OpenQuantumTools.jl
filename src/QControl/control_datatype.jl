abstract type DEDataControl <: AbstractAnnealingControl end
get_symbol(control::DEDataControl) = control.sym

require_de_data(::DEDataControl) = true
require_de_data(::AbstractAnnealingControl) = false
require_de_data(::Nothing) = false

function check_de_data_error(
    u0,
    control,
    de_data_constructor;
    additional_symbol = [],
)
    if require_de_data(control)
        if de_data_constructor == nothing
            error("Need to specify `de_array_constructor` for controller using `DEDataArray`.")
        else
            symbol_list = [:x]
            if !isempty(additional_symbol)
                symbol_list = vcat(symbol_list, additional_symbol)
            end
            controller_symbols = get_symbol(control)
            if isa(controller_symbols, Symbol)
                push!(symbol_list, controller_symbols)
            else
                symbol_list = vcat(symbol_list, controller_symbols)
            end
            for sym in symbol_list
                if !hasproperty(u0, sym)
                    error("No symbol :$controller_symbols defined in the DEDataArray.")
                end
            end
        end
    end
    nothing
end


(h::DenseHamiltonian)(du::DEDataArray, u::DEDataArray, tf, t) =
    (h::DenseHamiltonian)(du.x, u.x, tf, t)

QTBase.check_positivity(data::DEDataArray) = check_positivity(data.x)

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
