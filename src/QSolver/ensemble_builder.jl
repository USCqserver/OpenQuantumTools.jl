function build_ensemble_problem(
    A::Annealing,
    tf,
    type;
    span_unit = false,
    output_func = (sol, i) -> (sol, false),
    reduction = (u, data, I) -> (append!(u, data), false),
)
    if type == :stochastic_schrodinger
        prob, callback = build_ensemble_problem_stochastic_schrodinger(
            A,
            tf,
            DEFAULT_FLUCTUATOR_CONTROL_PROB_FUNC,
            output_func,
            reduction;
            span_unit = span_unit
        )
    else
        error("Ensemble problem of type $type is not implemented.")
    end
    prob, callback
end
