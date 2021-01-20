"""
$(SIGNATURES)

Build `EnsembleProblem` object for different open-system models. For `:lindblad` and `:ame`, it build the ensemble for the quantum trajectories method. For `:stochastic`, it builds the ensemble for the stochastic Schrödinger equation. For `:redfield`, it builds the ensemble to infuse stochastic nosie into the Redfield equation.

...
# Arguments
- `A::Annealing`: the `Annealing`/`Evolution` object.
- `tf::Real`: the total annealing time.
- `type::Symbol`: type of the ensemble to build. Available options are `:lindblad`, `:ame`, `:stochastic` and `:redfield`.
- `tspan = (0, tf)`: time interval to solve the dynamics.
- `initializer = DEFAULT_INITIALIZER`: initializer for the ensemble problem. Currently it is only supported by the `:stochastic` ensemble.
- `output_func`: The function determines what is saved from the solution to the output array. Defaults to saving the solution itself. The output is (out,rerun) where out is the output and rerun is a boolean which designates whether to rerun. It is part of the `DifferentialEquations`'s parallel interface.
- `reduction`: This function determines how to reduce the data in each batch. Defaults to appending the data from the batches. The second part of the output determines whether the simulation has converged. If true, the simulation will exit early. By default, this is always false. It is part of the `DifferentialEquations`'s parallel interface.
- `kwargs`: other keyword arguments supported by the specific solver or by `DifferentialEquations`.
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
    save_positions = (false, false),
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
            save_positions = save_positions,
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
            save_positions = save_positions,
            kwargs...,
        )
    elseif type == :lindblad
        res = build_ensemble_lindblad(
            A,
            tf,
            output_func,
            reduction;
            tspan = tspan,
            save_positions = save_positions,
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
    save_positions = (false, false),
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
            initializer = initializer,
            save_positions = save_positions,
            kwargs...,
        )
    else
        error("Ensemble problem of type $type is not implemented.")
    end
end

DEFAULT_PROB_FUNC(prob, i, repeat) = (prob)
