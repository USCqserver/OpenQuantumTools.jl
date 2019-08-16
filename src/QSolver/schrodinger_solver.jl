function solve_schrodinger(A::Annealing, tf::Real; span_unit = false, kwargs...)
    if ndims(A.u0) != 1
        throw(DomainError("Initial state must be given as a state vector."))
    end
    u0 = prepare_u0(A.u0, A.control)
    tf = prepare_tf(tf, span_unit)
    jp = sch_jacobian_prototype(A.H)
    p = AnnealingParams(A.H, tf; control = A.control)
    # tstops
    tstops = A.tstops
    if typeof(A.control) <: PausingControl
        ff = ODEFunction(
            sch_control_f;
            jac = sch_control_jac, jac_prototype = jp
        )
        cb = DiscreteCallback(pause_condition, pause_affect!)
        kwargs = Dict{Symbol,Any}(kwargs)
        kwargs[:callback] = cb
    else
        ff = ODEFunction(sch_f; jac = sch_jac, jac_prototype = jp)
    end
    tspan, tstops = scaling_time(tf, A.sspan, tstops)
    prob = ODEProblem{true}(ff, u0, tspan, p)
    solve(prob; alg_hints = [:nonstiff], tstops = tstops, kwargs...)
end

function solve_schrodinger(
    A::Annealing,
    tf::Vector{T},
    alg,
    para_alg = EnsembleSerial();
    output_func = (sol, i) -> (sol, false), span_unit = false, kwargs...
) where T <: Real
    if ndims(A.u0) != 1
        throw(DomainError("Initial state must be given as a state vector."))
    end
    u0 = prepare_u0(A.u0, A.control)
    t0 = prepare_tf(1.0, span_unit)
    jp = sch_jacobian_prototype(A.H)
    p = AnnealingParams(A.H, t0; control = A.control)
    # trajectories numbers
    trajectories = length(tf)
    tf_arr = float.(tf)
    # resolve control
    if typeof(A.control) <: PausingControl
        ff = ODEFunction(
            sch_control_f;
            jac = sch_control_jac, jac_prototype = jp
        )
        cb = DiscreteCallback(pause_condition, pause_affect!)
        kwargs = Dict{Symbol,Any}(kwargs)
        kwargs[:callback] = cb
    else
        ff = ODEFunction(sch_f; jac = sch_jac, jac_prototype = jp)
    end
    #
    if span_unit == true
        tstops = hyper_tstops(tf_arr, A.tstops)
        prob_func = (prob, i, repeat) -> begin
            tspan = (prob.tspan[1] * tf_arr[i], prob.tspan[2] * tf_arr[i])
            p = set_tf(prob.p, tf_arr[i])
            ODEProblem{true}(prob.f, prob.u0, tspan, p)
        end
    else
        tstops = A.tstops
        prob_func = (prob, i, repeat) -> begin
            p = set_tf(prob.p, tf_arr[i])
            ODEProblem{true}(prob.f, prob.u0, prob.tspan, p)
        end
    end
    prob = ODEProblem{true}(ff, u0, A.sspan, p)
    ensemble_prob = EnsembleProblem(
        prob;
        prob_func = prob_func, output_func = output_func
    )
    solve(
        ensemble_prob,
        alg,
        para_alg;
        trajectories = trajectories, tstops = tstops, kwargs...
    )
end


function sch_f(du, u, p::AbstractAnnealingParams, t::Real)
    p.H(du, u, p.tf, t)
end


function sch_jac(J, u, p, t)
    hmat = p.H(p.tf, t)
    mul!(J, -1.0im, hmat)
end


function sch_control_f(
    du::DEDataVector,
    u::DEDataVector,
    p::AbstractAnnealingParams,
    t::Real
)
    p.H(du, u, p, t)
end


function sch_control_jac(
    J,
    u::DEDataVector,
    p::AbstractAnnealingParams,
    t::Real
)
    hmat = p.H(u, p, t)
    mul!(J, -1.0im, hmat)
end
