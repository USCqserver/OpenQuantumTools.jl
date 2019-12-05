function build_ensemble_problem(
    A::Annealing,
    tf,
    type;
    span_unit = false,
    output_func = DEFAULT_OUTPUT_FUNC,
    prob_func = DEFAULT_PROB_FUNC,
    reduction = (u, data, I) -> (append!(u, data), false),
)
    if type == :stochastic_schrodinger
        if prob_func == DEFAULT_PROB_FUNC
            prob_func = DEFAULT_FLUCTUATOR_CONTROL_PROB_FUNC
        end
        prob, callback = build_ensemble_problem_stochastic_schrodinger(
            A,
            tf,
            prob_func,
            output_func,
            reduction;
            span_unit = span_unit,
        )
    elseif type == :schrodinger
        if prob_func == DEFAULT_PROB_FUNC
            tf_arr = ndims(tf) == 0 ?
                     error("Default schrodinger ensemble model parallelized for tf array.") :
                     float.(tf)
            prob_func = (prob, i, repeat) -> begin
                p = set_tf(prob.p, tf_arr[i])
                ODEProblem{true}(prob.f, prob.u0, prob.tspan, p)
            end
        end
        prob, callback = build_ensemble_problem_schrodinger(
            A,
            prob_func,
            output_func,
            reduction;
            span_unit = span_unit,
        )
    else
        error("Ensemble problem of type $type is not implemented.")
    end
    prob, callback
end

DEFAULT_OUTPUT_FUNC(sol, i) = (sol, false)
# this is not the actually default problem function
DEFAULT_PROB_FUNC(prob, i, repeat) = (prob)
