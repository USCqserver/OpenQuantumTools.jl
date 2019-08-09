function solve_schrodinger(A::Annealing, tf::Real; kwargs...)
    if ndims(A.u0) != 1
        throw(DomainError("Initial state must be given as a state vector."))
    end
    u0 = sch_prepare_u0(A.u0)
    p = AnnealingParams(A.H, float(tf))
    jp = zeros(eltype(u0), length(u0), length(u0))
    if haskey(kwargs, :span_unit) && kwargs[:span_unit] == true
        t_span = (tf * A.sspan[1], tf * A.sspan[2])
        ff = ODEFunction(sch_f_unit; jac = sch_jac_unit, jac_prototype = jp)
    else
        t_span = A.sspan
        ff = ODEFunction(
            sch_f_unitless;
            jac = sch_jac_unitless, jac_prototype = jp
        )
    end
    prob = ODEProblem{true}(ff, u0, t_span, p)
    solve(prob; alg_hints = [:nonstiff], tstops = A.tstops, kwargs...)
end

function solve_schrodinger(
    A::Annealing,
    tf::Vector{T},
    alg,
    para_alg = EnsembleSerial();
    kwargs...
) where T <: Real
    if ndims(A.u0) != 1
        throw(DomainError("Initial state must be given as a state vector."))
    end
    u0 = sch_prepare_u0(A.u0)
    p = AnnealingParams(A.H, float(tf[1]))
    jp = zeros(eltype(u0), length(u0), length(u0))
    # trajectories numbers
    trajectories = length(tf)
    tf_arr = float.(tf)
    if haskey(kwargs, :span_unit) && kwargs[:span_unit] == true
        t_span = (float(tf[1]) * A.sspan[1], float(tf[1]) * A.sspan[2])
        ff = ODEFunction(sch_f_unit; jac = sch_jac_unit, jac_prototype = jp)
    else
        t_span = A.sspan
        ff = ODEFunction(
            sch_f_unitless;
            jac = sch_jac_unitless, jac_prototype = jp
        )
    end
    function prob_func(prob, i, repeat)
        prob.p.H = p_copy(prob.p.H)
        prob.p.tf = tf_arr[i]
        prob
    end
    prob = ODEProblem{true}(ff, u0, t_span, p)
    ensemble_prob = EnsembleProblem(prob; prob_func = prob_func)
    solve(
        ensemble_prob,
        alg,
        para_alg;
        trajectories = trajectories, tstops = A.tstops, kwargs...
    )
end

function sch_f_unitless(du, u, p, t)
    p.H(du, u, p.tf, t)
end

function sch_f_unit(du, u, p, t)
    p.H(du, u, 1.0, t / p.tf)
end

function sch_jac_unitless(J, u, p, t)
    hmat = p.H(t)
    mul!(J, -1.0im * p.tf, hmat)
end

function sch_jac_unit(J, u, p, t)
    hmat = p.H(t / p.tf)
    mul!(J, -1.0im, hmat)
end

function sch_prepare_u0(raw_u0, control = false)
    if control == false
        complex(raw_u0)
    else
        SchVec(complex(raw_u0), 1.0)
    end
end

mutable struct SchVec{T} <: DEDataVector{T}
    x::Array{T,1}
    state::T
end
