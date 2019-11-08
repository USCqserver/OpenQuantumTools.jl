"""
    function prepare_u0(raw_u0; type =:v, control=nothing, vectorize=false)

Prepare initial state in proper type and dimension for ODE solvers. `type` specifies the dimension of the initial state: `:v` is 1-D state vector and `:m` is 2-D density matrix. `control` should be any `AbstractAnnealingControl` object. `vectorize` indicate whether to vectorize the density matrix.
"""
function prepare_u0(raw_u0; type =:v, control=nothing, vectorize=false)
    res = complex(raw_u0)
    if type == :v && ndims(res) != 1
        throw(ArgumentError("Cannot convert density matrix to state vector."))
    elseif type ==:m && ndims(res) == 1
        res = res*res'
    elseif ndims(res) < 1 || ndims(res) > 2
        throw(ArgumentError("u0 can either be a vector or matrix."))
    end
    if control != nothing
        res = adjust_u0_with_control(res, control)
    end
    if vectorize == true
        res = res[:]
    end
    res
end


function create_tstops_for_tf_array(tf_arr, tstops)
    temp = [tf * tstops for tf in tf_arr]
    temp = vcat(temp...)
    unique(sort(temp))
end


function prepare_tf(tf, span_unit)
    if span_unit == true
        UnitTime(tf)
    else
        float(tf)
    end
end


function scaling_time(tf::UnitTime, tspan, tstops)
    (tf * tspan[1], tf * tspan[2]), tf * tstops
end


function scaling_time(tf::Real, tspan, tstops)
    tspan, tstops
end


function vectorized_jacobian_prototype(H)
    num_type = typeof(H).parameters[1]
    iden = Matrix{num_type}(I, H.size)
    if typeof(H) <: AdiabaticFrameHamiltonian
        iden ⊗ zeros(num_type, H.size)
    else
        iden ⊗ H.u_cache
    end
end


function sch_jacobian_prototype(H)
    if typeof(H) <: AdiabaticFrameHamiltonian
        zeros(typeof(H).parameters[1], H.size)
    else
        similar(H.u_cache)
    end
end


function prepare_callback(kw_dict, control)
    kw_dict
end
