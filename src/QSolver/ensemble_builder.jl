"""
$(SIGNATURES)

Build `EnsembleProblem` object for different opensystem dynamics of `Annealing` process `A`.

...
# Arguments
- `A::Annealing`: the Annealing object.
- `tf::Real`: the total annealing time.
- `type::Symbol`: type of master equation to build. Available options are `:ame`, `:stochastic` and `:redfield`.
- `tspan = (0, tf)`: time interval to solve.
- `output_func`
- `reduction`
- `kwargs`: other keyword arguments supported by DifferentialEquations.jl.
...
"""
function build_ensembles(
    A::Annealing,
    tf::Real,
    type::Symbol;
    tspan = (0.0, tf),
    output_func = (sol, i) -> (sol, false),
    reduction = (u, data, I) -> (append!(u, data), false),
    initializer = DEFAULT_INITIALIZER,
    kwargs...,
)
    if type == :stochastic
        res = build_ensemble_stochastic(
            A,
            tf,
            output_func,
            reduction;
            tspan = tspan,
            initializer = initializer,
            kwargs...,
        )
    elseif type == :ame
        res = build_ensemble_ame(
            A,
            tf,
            output_func,
            reduction;
            tspan = tspan,
            initializer = initializer,
            kwargs...,
        )
    elseif type == :lindblad
        res = build_ensemble_lindblad(
            A,
            tf,
            output_func,
            reduction;
            tspan = tspan,
            kwargs...,
        )
    else
        error("Ensemble problem of type $type is not supported.")
    end
    res
end

function build_ensembles(
    A::Annealing,
    tf::Real,
    unitary,
    type::Symbol;
    vectorize::Bool = false,
    output_func = (sol, i) -> (sol, false),
    reduction = (u, data, I) -> (append!(u, data), false),
    int_atol = 1e-8,
    int_rtol = 1e-6,
    Ta = tf,
    initializer = DEFAULT_INITIALIZER,
    kwargs...,
)
    if type == :redfield
        build_ensemble_redfield(
            A,
            tf,
            unitary,
            output_func,
            reduction;
            vectorize = vectorize,
            int_atol = int_atol,
            int_rtol = int_rtol,
            Ta = Ta,
            initializer = DEFAULT_INITIALIZER,
            kwargs...,
        )
    else
        error("Ensemble problem of type $type is not implemented.")
    end
end

DEFAULT_PROB_FUNC(prob, i, repeat) = (prob)
