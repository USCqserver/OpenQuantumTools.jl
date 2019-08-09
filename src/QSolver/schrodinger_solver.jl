function solve_schrodinger(A::Annealing, tf::Real; kwargs...)
    if ndims(A.u0) != 1
        throw(DomainError("Initial state must be given as a state vector."))
    end
    u0 = sch_prepare_u0(A.u0)
    p = AnnealingParams(A.H, float(tf))
    jp = zeros(eltype(u0), length(u0), length(u0))
    if haskey(kwargs, :span_unit) && kwargs[:span_unit] == true
        tspan = (tf * A.sspan[1], tf * A.sspan[2])
        tstops = tf * A.tstops
        ff = ODEFunction(sch_f_unit; jac = sch_jac_unit, jac_prototype = jp)
    else
        tspan = A.sspan
        tstops = A.tstops
        ff = ODEFunction(
            sch_f_unitless;
            jac = sch_jac_unitless, jac_prototype = jp
        )
    end
    prob = ODEProblem{true}(ff, u0, tspan, p)
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
    p = AnnealingParams(A.H, 1.0)
    jp = zeros(eltype(u0), length(u0), length(u0))
    # trajectories numbers
    trajectories = length(tf)
    tf_arr = float.(tf)
    if haskey(kwargs, :span_unit) && kwargs[:span_unit] == true
        if !isempty(A.tstops)
            @warn "Parallel algorithm does not support scaling tstops argument. Setting it to be empty."
            tstops = []
        end
        ff = ODEFunction(sch_f_unit; jac = sch_jac_unit, jac_prototype = jp)
        prob_func = (
            prob,
            i,
            repeat
        ) -> begin
            prob.p.H = p_copy(prob.p.H)
            tspan = (
                prob.tspan[1] * tf_arr[i],
                prob.tspan[2] * tf_arr[i]
            )
            prob.p.tf = tf_arr[i]
            ODEProblem{true}(prob.f, prob.u0, tspan, prob.p)
        end
    else
        tstops = A.tstops
        ff = ODEFunction(
            sch_f_unitless;
            jac = sch_jac_unitless, jac_prototype = jp
        )
        prob_func = (
            prob,
            i,
            repeat
        ) -> begin
            prob.p.H = p_copy(prob.p.H)
            prob.p.tf = tf_arr[i]
            prob
        end
    end
    prob = ODEProblem{true}(ff, u0, A.sspan, p)
    ensemble_prob = EnsembleProblem(prob; prob_func = prob_func)
    solve(
        ensemble_prob,
        alg,
        para_alg;
        trajectories = trajectories, tstops = A.tstops, kwargs...
    )
end

function ode_with_unit() end

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
