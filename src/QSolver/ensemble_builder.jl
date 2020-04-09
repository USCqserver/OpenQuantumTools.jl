function build_ensemble_problem(
    A::Annealing,
    tf,
    type;
    dimensionless_time = true,
    output_func = DEFAULT_OUTPUT_FUNC,
    reduction = (u, data, I) -> (append!(u, data), false),
    de_array_constructor = nothing,
    additional_field = nothing,
    initializer =  DEFAULT_INITIALIZER,
    kwargs...
)
    if type == :stochastic_schrodinger
        res = build_ensemble_problem_stochastic_schrodinger(
            A,
            tf,
            output_func,
            reduction;
            dimensionless_time = dimensionless_time,
            de_array_constructor = de_array_constructor,
            fluctuator_de_field = additional_field,
            initializer = initializer
        )
    elseif type == :ame_trajectory
        if prob_func == DEFAULT_PROB_FUNC
            prob_func = DEFAULT_AME_TRAJECTORY_PROB_FUNC
        end
    else
        error("Ensemble problem of type $type is not implemented.")
    end
    # builder = getfield(QuantumAnnealingTools, Symbol("build_ensemble_problem_", type))
    # prob, callback = builder(
    #     A,
    #     tf,
    #     prob_func,
    #     output_func,
    #     reduction;
    #     dimensionless_time = dimensionless_time,
    #     de_array_constructor = de_array_constructor,
    #     fluctuator_de_field = fluctuator_de_field,
    #     kwargs...
    # )
    res
end

DEFAULT_OUTPUT_FUNC(sol, i) = (sol, false)
# this is not the actually default problem function
DEFAULT_PROB_FUNC(prob, i, repeat) = (prob)


function build_prob_func(func, type::Symbol)
    if type == :fluctuator
        res = function (prob, i, repeat)
            ctrl = prob.p.control
            ctrl.b0 = abs.(ctrl.b0) .* func()
            u0 = prob.u0
            u0.n .= ctrl()
            next_state!(ctrl)
            ODEProblem{true}(prob.f, u0, prob.tspan, prob.p)
        end
    else
        error("Build prob_func for type $type is not implemented.")
    end
    res
end


function DEFAULT_AME_TRAJECTORY_PROB_FUNC(prob, i, repeat)
    prob.u0.r = rand()
    prob
end
