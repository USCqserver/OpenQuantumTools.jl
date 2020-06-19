"""
This function initialize `(tf, u0, tstops)` based on all the other input arguments of a solver.
"""
function __init(u0, tf, dimensionless, s_type, tstops, a_tstops, control,
    de_builder; vectorize = false, needed_symbol = [])
    #TODO: Refactor this to a better format
    # If tf is 0, always use dimensionless time
    # This is to keep the solution format uniform and avoid t/tf
    if tf == 0
        dimensionless = true
    end
    u0 = build_u0(u0, s_type, vectorize, de_builder)
    tf, tstops = preprocessing_time(tf, tstops, a_tstops, dimensionless)
    check_de_data_error(u0, control, de_builder, needed_symbol)
    tf, u0, tstops
end

__init(A::Annealing, tf, dimensionless, s_type, tstops, de_builder;
    vectorize = false, needed_symbol = []) = __init(A.u0, tf, dimensionless,
    s_type, tstops, A.tstops, A.control, de_builder; vectorize = vectorize,
    needed_symbol = needed_symbol)

"""
    function preprocessing_time(tf, tstops_in_args, tstops_pre_defined, dimensionless)

Preprocessing the total evolution time `tf` and the extra time stepping points `tstops` depending on `dimensionless` argument. The function will combine `tstops_in_args` and `tstops_pre_defined`.
"""
function preprocessing_time(
    tf,
    tstops_in_args,
    tstops_pre_defined,
    dimensionless,
)
    tf = dimensionless == true ? float(tf) : UnitTime(tf)
    if !isempty(tstops_in_args) && !isempty(tstops_pre_defined)
        @warn "Both the annealing object and the solver arguments have tstops. They will be merged together."
    end
    tstops = vcat(tstops_in_args, tstops_pre_defined)
    tstops = typeof(tf) <: UnitTime ? tf * tstops : tstops
    tf, tstops
end

"""
    function build_u0(raw_u0, type, vectorize, de_builder)

Prepare initial state in proper type and dimension for ODE solvers. `type` specifies the dimension of the initial state: `:v` is 1-D state vector and `:m` is 2-D density matrix. `vectorize` indicate whether to vectorize the density matrix. `de_builder` specify whether to convert the resulting `u0` into a `DEDataArray` type.
"""
function build_u0(raw_u0, type, vectorize, de_builder)
    res = complex(raw_u0)
    if type == :v && ndims(res) != 1
        throw(ArgumentError("Cannot convert density matrix to state vector."))
    elseif type == :m && ndims(res) == 1
        res = res * res'
    elseif ndims(res) < 1 || ndims(res) > 2
        throw(ArgumentError("u0 can either be a vector or matrix."))
    end
    if vectorize == true
        res = res[:]
    end
    de_builder == nothing ? res : de_builder(res)
end

function check_de_data_error(u0, control, de_builder, additional_symbol)
    if require_de_data(control)
        if de_builder == nothing
            error("Need to specify `de_array_constructor` for controller using `DEDataArray`.")
        else
            symbol_list = [:x]
            if !isempty(additional_symbol)
                symbol_list = vcat(symbol_list, additional_symbol)
            end
            controller_symbols = get_symbol(control)
            if isa(controller_symbols, Symbol)
                push!(symbol_list, controller_symbols)
            else
                symbol_list = vcat(symbol_list, controller_symbols)
            end
            for sym in symbol_list
                if !hasproperty(u0, sym)
                    error("No symbol :$controller_symbols defined in the DEDataArray.")
                end
            end
        end
    end
    nothing
end

scaling_tspan(tf::Real, tspan) = tspan
scaling_tspan(tf::UnitTime, tspan) = (tf * tspan[1], tf * tspan[2])

function build_prob_func(initializer)
    prob_func = function (prob, i, repeat)
        ctrl = prob.p.control
        reset!(ctrl, prob.u0, initializer)
        ODEProblem{true}(prob.f, prob.u0, prob.tspan, prob.p)
    end
end
