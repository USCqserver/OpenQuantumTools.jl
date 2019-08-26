function prepare_u0(raw_u0, control)
    res = complex(raw_u0)
    if typeof(control) <: PausingControl
        if ndims(raw_u0) == 1
            res = DEPausingVec(raw_u0, 1)
        elseif ndims(raw_u0) == 2
            res = DEPausingMat(raw_u0, 1)
        else
            throw(ArgumentError("u0 can either be a vector or matrix."))
        end
    end
    res
end


function hyper_tstops(tf_arr, tstops)
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
