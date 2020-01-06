function build_ensemble_problem(
    A::Annealing,
    tf,
    type;
    span_unit = false,
    output_func = DEFAULT_OUTPUT_FUNC,
    prob_func = DEFAULT_PROB_FUNC,
    reduction = (u, data, I) -> (append!(u, data), false),
    kwargs...
)
    if type == :stochastic_schrodinger
        if prob_func == DEFAULT_PROB_FUNC
            prob_func = DEFAULT_FLUCTUATOR_CONTROL_PROB_FUNC
        end
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
    elseif type == :ame_trajectory
        if prob_func == DEFAULT_PROB_FUNC
            prob_func = DEFAULT_AME_TRAJECTORY_PROB_FUNC
        end
    else
        error("Ensemble problem of type $type is not implemented.")
    end
    builder = getfield(QuantumAnnealingTools, Symbol("build_ensemble_problem_", type))
    prob, callback = builder(
        A,
        tf,
        prob_func,
        output_func,
        reduction;
        span_unit = span_unit,
        kwargs...
    )
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
