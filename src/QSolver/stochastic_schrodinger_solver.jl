function solve_stochastic_schrodinger(
    A::Annealing,
    tf::Real;
    dimensionless_time = true,
    tstops = Float64[],
    de_array_constructor = nothing,
    fluctuator_de_field = nothing,
    initializer = DEFAULT_INITIALIZER,
    kwargs...,
)
    tf, tstops = preprocessing_time(tf, tstops, A.tstops, dimensionless_time)
    u0 = build_u0(A.u0, :v, de_array_constructor = de_array_constructor)
    additional_symbol =
        fluctuator_de_field == nothing ? [] : [fluctuator_de_field]
    check_de_data_error(
        u0,
        A.control,
        de_array_constructor,
        additional_symbol = additional_symbol,
    )
    # build control object from bath; a prototype implementation
    if A.interactions == nothing
        control = FluctuatorControl(
            tf,
            length(A.coupling),
            A.bath,
            fluctuator_de_field,
        )
        opensys = StochasticNoise(A.coupling, fluctuator_de_field)
    else
        error("Interactions are not supported for stochastic schrodinger equation yet.")
    end
    reset!(A.control)
    reset!(control, u0, initializer)
    control = A.control == nothing ? control : ControlSet(control, A.control)

    callback = stochastic_schrodinger_build_callback(control)
    p = ODEParams(A.H, tf; opensys = opensys, control = control)
    ff = stochastic_schrodinger_build_ode_function(A.H, control)
    prob = ODEProblem{true}(ff, u0, (p) -> scaling_tspan(p.tf, A.sspan), p)

    solve(
        prob;
        alg_hints = [:nonstiff],
        callback = callback,
        tstops = tstops,
        kwargs...,
    )
end


function build_ensemble_problem_stochastic_schrodinger(
    A::Annealing,
    tf::Real,
    output_func,
    reduction;
    dimensionless_time = true,
    de_array_constructor = nothing,
    fluctuator_de_field = nothing,
    initializer = DEFAULT_INITIALIZER,
    tstops = Float64[],
)
    tf, tstops = preprocessing_time(tf, tstops, A.tstops, dimensionless_time)
    u0 = build_u0(A.u0, :v, de_array_constructor = de_array_constructor)
    additional_symbol =
        fluctuator_de_field == nothing ? [] : [fluctuator_de_field]
    check_de_data_error(
        u0,
        A.control,
        de_array_constructor,
        additional_symbol = additional_symbol,
    )
    # build control object from bath; a prototype implementation
    if A.interactions == nothing
        control = FluctuatorControl(
            tf,
            length(A.coupling),
            A.bath,
            fluctuator_de_field,
        )
        opensys = StochasticNoise(A.coupling, fluctuator_de_field)
    else
        error("Interactions are not supported for stochastic schrodinger equation yet.")
    end

    control = A.control == nothing ? control : ControlSet(control, A.control)
    callback = stochastic_schrodinger_build_callback(control)
    prob_func = build_stochastic_prob_func(control, initializer)
    p = ODEParams(A.H, tf; opensys = opensys, control = control)
    ff = stochastic_schrodinger_build_ode_function(A.H, control)
    prob = ODEProblem{true}(ff, u0, (p) -> scaling_tspan(p.tf, A.sspan), p)

    ensemble_prob = EnsembleProblem(
        prob;
        prob_func = prob_func,
        output_func = output_func,
        reduction = reduction,
    )

    ensemble_prob, callback
end


function stochastic_schrodinger_build_callback(control::ControlSet)
    callbacks = [stochastic_schrodinger_build_callback(v, k) for (k, v) in zip(keys(control), control)]
    CallbackSet(callbacks...)
end


stochastic_schrodinger_build_callback(control::Union{FluctuatorControl, FluctuatorDEControl}) =
    build_callback(control)


stochastic_schrodinger_build_callback(control::Union{FluctuatorControl, FluctuatorDEControl}, sym::Symbol) =
    build_callback(control, sym)


stochastic_schrodinger_build_callback(
    control::Union{InstPulseControl,InstDEPulseControl},
    sym::Symbol,
) = build_callback(control, sym, (c, pulse) -> c .= pulse * c)


function stochastic_schrodinger_build_ode_function(H, control)
    update_func = build_update_func(control)
    cache = get_cache(H)
    diff_op = DiffEqArrayOperator(cache, update_func = update_func)
    jac_cache = similar(cache)
    jac_op = DiffEqArrayOperator(jac_cache, update_func = update_func)
    ff = ODEFunction(diff_op; jac_prototype = jac_op)
end

function build_update_func(::FluctuatorControl)
    update_func = function (A, u, p, t)
        update_cache!(A, p.H, p.tf, t)
        p.opensys(A, p.control(), p.tf, t)
    end
end


function build_update_func(::FluctuatorDEControl)
    update_func = function (A, u, p, t)
        update_cache!(A, p.H, p.tf, t)
        p.opensys(A, u, p.tf, t)
    end
end


function build_update_func(control::ControlSet)
    f_control = control.fluctuator_control
    if typeof(f_control) <: FluctuatorControl
        update_func = function (A, u, p, t)
            update_cache!(A, p.H, p.tf, t)
            c = p.control.fluctuator_control
            p.opensys(A, c.n, p.tf, t)
        end
    else
        update_func = function (A, u, p, t)
            update_cache!(A, p.H, p.tf, t)
            p.opensys(A, u, p.tf, t)
        end
    end
    update_func
end


function build_stochastic_prob_func(::Union{FluctuatorControl, FluctuatorDEControl}, initializer)
    prob_func = function (prob, i, repeat)
        ctrl = prob.p.control
        reset!(ctrl, prob.u0, initializer)
        ODEProblem{true}(prob.f, prob.u0, prob.tspan, prob.p)
    end
end


function build_stochastic_prob_func(::ControlSet, initializer)
    prob_func = function (prob, i, repeat)
        ctrl = getproperty(prob.p.control, :fluctuator_control)
        reset!(ctrl, prob.u0, initializer)
        ODEProblem{true}(prob.f, prob.u0, prob.tspan, prob.p)
    end
end
