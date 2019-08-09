function solve_davies(A::Annealing, tf::Real; kwargs...)
    H = A.H
    # turn initial state to density matrix
    if ndims(A.u0) == 1
        u0 = A.u0*A.u0'
    else
        u0 = A.u0
    end
    # check if an estimation of frequency range is given
    if haskey(kwargs, :ω_hint)
        ωr = kwargs[:ω_hint]
    else
        ωr = nothing
    end
    #
    if haskey(kwargs, :lvl)
        lvl = kwargs[:lvl]
    else
        lvl = A.H.size[1]
    end
    #
    if typeof(A.H) <: AbstractSparseHamiltonian
        if A.H.size[1] == 2
            @warn "Hamiltonian size is too small for sparse factorization. Convert to dense Hamiltonian"
            H = to_dense(A.H)
        elseif lvl == A.H.size[1]
            @warn "Sparse Hamiltonian detected. Truncate the level to n-1."
            lvl = lvl-1
        end
    end
    γ, S = davies_spectrum(A.bath; ω_range = ωr)
    DiffOp = DaviesDiffEqOperator()
    p = AnnealingParams(A.H, float(tf); opensys=opensys)
    prob = ODEProblem(davies_ode, u0, A.sspan, p)
    solve(prob; alg_hints = [:nonstiff], tstops=A.tstops, kwargs...)
end


function solve_davies(A::Annealing, tf::Vector{T}, para_alg = EnsembleSerial(); kwargs...) where T<:Number
    if ndims(A.u0) == 1
        u0 = A.u0*A.u0'
    else
        u0 = A.u0
    end
    if haskey(kwargs, :ω_hint)
        ωr = kwargs[:ω_hint]
    else
        ωr = nothing
    end
    opensys = create_davies(A.coupling, A.bath; ω_range = ωr)
    trajectories = length(tf)
    tf_arr = float.(tf)
    p = AnnealingParams(A.H, float(tf[1]); opensys=opensys)
    prob = ODEProblem(davies_ode, u0, A.sspan, p)
    function prob_func(prob, i, repeat)
        prob.p.tf = tf_arr[i]
        prob.p.H = p_copy(prob.p.H)
        prob
    end
    ensemble_prob = EnsembleProblem(prob, prob_func=prob_func)
    solve(ensemble_prob, para_alg; alg_hints = [:nonstiff], tstops=A.tstops, trajectories=trajectories, kwargs...)
end


function davies_spectrum(bath::OhmicBath; ω_range)
    if ω_range == nothing
        γ_loc(ω) = γ(ω, bath)
        S_loc(ω) = S(ω, bath)
    else
        γ_loc, S_loc = interpolate_spectral_density(2π*ω_range, bath)
    end
    γ_loc, S_loc
end

function create_davies(coupling, bath::OhmicBath; ω_range = nothing)
    if ω_range == nothing
        γ_loc(ω) = γ(ω, bath)
        S_loc(ω) = S(ω, bath)
    else
        γ_loc, S_loc = interpolate_spectral_density(ω_range, bath)
    end
    Davies(coupling, γ_loc, S_loc)
end
