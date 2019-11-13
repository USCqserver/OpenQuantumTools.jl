"""
    function prepare_u0(raw_u0; type =:v, control=nothing, vectorize=false)

Prepare initial state in proper type and dimension for ODE solvers. `type` specifies the dimension of the initial state: `:v` is 1-D state vector and `:m` is 2-D density matrix. `control` should be any `AbstractAnnealingControl` object. `vectorize` indicate whether to vectorize the density matrix.
"""
function prepare_u0(raw_u0; type = :v, control = nothing, vectorize = false)
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
    if control != nothing
        res = adjust_u0_with_control(res, control)
    end
    res
end


function prepare_tspan(tspan)
    (p)->scaling_tspan(p.tf, tspan)
end


function prepare_tstops(tf, arg_tstops, obj_tstops)
    if !isempty(arg_tstops) && !isempty(obj_tstops)
        @warn "Both the annealing object and the solver arguments have tstops. They will be merged together."
    end
    tstops = vcat(arg_tstops, obj_tstops)
    if typeof(tf) <: UnitTime
        return tf * tstops
    else
        return tstops
    end
end


function create_tstops_for_tf_array(tf_arr, tstops)
    temp = [tf * tstops for tf in tf_arr]
    temp = vcat(temp...)
    unique(sort(temp))
end


prepare_tf(tf, span_unit) = span_unit==true ? UnitTime(tf) : float(tf)
scaling_tspan(tf::Real, tspan) = tspan
scaling_tspan(tf::UnitTime, tspan) = (tf*tspan[1], tf*tspan[2])

scaling_time(tf::UnitTime, tspan, tstops) = (tf*tspan[1], tf*tspan[2]), tf*tstops
scaling_time(tf::Real, tspan, tstops) = tspan, tstops
