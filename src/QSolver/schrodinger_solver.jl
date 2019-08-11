function solve_schrodinger(A::Annealing, tf::Real; kwargs...)
    if ndims(A.u0) != 1
        throw(DomainError("Initial state must be given as a state vector."))
    end
    u0 = prepare_u0(A.u0)
    jp = zeros(eltype(u0), length(u0), length(u0))
    if haskey(kwargs, :span_unit) && kwargs[:span_unit] == true
        tspan = (tf * A.sspan[1], tf * A.sspan[2])
        tstops = tf * A.tstops
        param_tf = UnitTime(tf)
    else
        tspan = A.sspan
        tstops = A.tstops
        param_tf = float(tf)
    end
    p = AnnealingParams(A.H, param_tf)
    ff = ODEFunction(sch_f; jac = sch_jac, jac_prototype = jp)
    prob = ODEProblem{true}(ff, u0, tspan, p)
    solve(prob; alg_hints = [:nonstiff], tstops = tstops, kwargs...)
end

function solve_schrodinger(
    A::Annealing,
    tf::Vector{T},
    alg,
    para_alg = EnsembleSerial();
    output_func = (sol, i) -> (sol, false), kwargs...
) where T <: Real
    if ndims(A.u0) != 1
        throw(DomainError("Initial state must be given as a state vector."))
    end
    u0 = prepare_u0(A.u0)
    jp = zeros(eltype(u0), length(u0), length(u0))
    # trajectories numbers
    trajectories = length(tf)
    tf_arr = float.(tf)
    if haskey(kwargs, :span_unit) && kwargs[:span_unit] == true
        if !isempty(A.tstops)
            @warn "Parallel algorithm does not support scaling tstops argument. Setting it to be empty."
            tstops = []
        end
        p = AnnealingParams(A.H, UnitTime(1.0))
        prob_func = (
            prob,
            i,
            repeat
        ) -> begin
            tspan = (prob.tspan[1] * tf_arr[i], prob.tspan[2] * tf_arr[i])
            p = set_tf(prob.p, tf_arr[i])
            ODEProblem{true}(prob.f, prob.u0, tspan, p)
        end
    else
        tstops = A.tstops
        p = AnnealingParams(A.H, 1.0)
        prob_func = (
            prob,
            i,
            repeat
        ) -> begin
            p = set_tf(prob.p, tf_arr[i])
            ODEProblem{true}(prob.f, prob.u0, prob.tspan, p)
        end
    end
    ff = ODEFunction(sch_f; jac = sch_jac, jac_prototype = jp)
    prob = ODEProblem{true}(ff, u0, A.sspan, p)
    ensemble_prob = EnsembleProblem(
        prob;
        prob_func = prob_func, output_func = output_func
    )
    solve(
        ensemble_prob,
        alg,
        para_alg;
        trajectories = trajectories, tstops = A.tstops, kwargs...
    )
end


function sch_f(du, u, p, t)
    p.H(du, u, p.tf, t)
end


function sch_jac(J, u, p, t)
    hmat = p.H(p.tf, t)
    mul!(J, -1.0im, hmat)
end
