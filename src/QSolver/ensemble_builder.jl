function build_ensemble_problem(
    A::Annealing,
    tf,
    type;
    dimensionless_time = true,
    output_func = DEFAULT_OUTPUT_FUNC,
    reduction = (u, data, I) -> (append!(u, data), false),
    de_array_constructor = nothing,
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
            initializer = initializer,
            kwargs...
        )
    elseif type == :ame_trajectory
        res = build_ensemble_problem_ame_trajectory(
            A,
            tf,
            output_func,
            reduction;
            dimensionless_time = dimensionless_time,
            de_array_constructor = de_array_constructor,
            initializer = initializer,
            kwargs...
        )
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
