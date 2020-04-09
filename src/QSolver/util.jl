"""
    function preprocessing_time(tf, tstops_in_args, tstops_pre_defined, dimensionless_time)

Preprocessing the total evolution time `tf` and the extra time stepping points `tstops` depending on `dimensionless_time` argument. The function will combine `tstops_in_args` and `tstops_pre_defined`.
"""
function preprocessing_time(
    tf,
    tstops_in_args,
    tstops_pre_defined,
    dimensionless_time,
)
    tf = dimensionless_time == true ? float(tf) : UnitTime(tf)
    if !isempty(tstops_in_args) && !isempty(tstops_pre_defined)
        @warn "Both the annealing object and the solver arguments have tstops. They will be merged together."
    end
    tstops = vcat(tstops_in_args, tstops_pre_defined)
    tstops = typeof(tf) <: UnitTime ? tf * tstops : tstops
    tf, tstops
end


"""
    function build_u0(raw_u0, type; vectorize=false)

Prepare initial state in proper type and dimension for ODE solvers. `type` specifies the dimension of the initial state: `:v` is 1-D state vector and `:m` is 2-D density matrix. `vectorize` indicate whether to vectorize the density matrix.
"""
function build_u0(
    raw_u0,
    type;
    vectorize = false,
    de_array_constructor = nothing,
)
    res = complex(raw_u0)
    if type == :v && ndims(res) != 1
        throw(ArgumentError("Cannot convert density matrix to state vector."))
    elseif type == :m && ndims(res) == 1
        res = res * res'
    elseif ndims(res) < 1 || ndims(res) > 2
        throw(ArgumentError("u0 can either be a vector or matrix."))
    end
    if vectorize == true
        res = res[:]
    end
    de_array_constructor == nothing ? res : de_array_constructor(res)
end

reset!(::Nothing) = nothing
scaling_tspan(tf::Real, tspan) = tspan
scaling_tspan(tf::UnitTime, tspan) = (tf * tspan[1], tf * tspan[2])


"""
$(TYPEDEF)
Defines a complete set of ODE parameters, which includes Hamiltonian, total annealing time, open system and control objects.
# Fields
$(FIELDS)
"""
struct ODEParams{T<:Union{AbstractFloat,UnitTime}}
    """Hamiltonian"""
    H
    """Total annealing time"""
    tf::T
    """Open system object"""
    opensys
    """Annealing control object"""
    control
end


function ODEParams(
    H,
    tf::T;
    opensys = nothing,
    control = nothing,
) where {T<:Number}
    ODEParams(H, float(tf), opensys, control)
end

ODEParams(H, tf::UnitTime; opensys = nothing, control = nothing) =
    ODEParams(H, tf, opensys, control)

ODEParams(tf; opensys = nothing, control = nothing) =
    ODEParams(nothing, tf, opensys, control)
