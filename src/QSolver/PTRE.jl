"""
    function SA_Δ²(t_idx, tf, TG, C, bath::HybridOhmicBath, i = 1, j =2)

Calculate the tunneling matrix elements of the Smirnov-Amin ME -- ``Δ^2_{ij}``.
"""
function SA_Δ²(t_idx, tf, TG, C, bath::HybridOhmicBath, i = 1, j = 2)
    # prepare parameters
    # TODO refractor out the following code block
    ω = 2π * (TG.ω[t_idx, j] - TG.ω[t_idx, i])
    a = C.a[t_idx, j, i]
    b = C.b[t_idx, j, i]
    c = C.c[t_idx, j, i]
    d = C.d[t_idx, j, i]
    T = TG.T[t_idx, j, i]
    G = TG.G[t_idx, j, i]

    T̃ = T - 1.0im * G / tf - d * bath.ϵ
    A = abs2(T̃ - ω * c / a)
    B = (a * b - abs2(c)) / a^2
    A + B * (ω^2 + a * bath.W^2)
end


"""
    function SA_redfield(t_idx, tf, TG, C, bath::HybridOhmicBath, i = 1, j =2)

Calculate the Redfield rate of Smirnov-Amin ME -- ``Γ_{ij}``.
"""
function SA_redfield(t_idx, tf, TG, C, bath::HybridOhmicBath, i = 1, j = 2)
    ω = 2π * (TG.ω[t_idx, j] - TG.ω[t_idx, i])
    a = C.a[t_idx, j, i]
    b = C.b[t_idx, j, i]
    c = C.c[t_idx, j, i]
    d = C.d[t_idx, j, i]
    T = TG.T[t_idx, j, i]
    G = TG.G[t_idx, j, i]

    T̃ = T - 1.0im * G / tf - d * bath.ϵ
    A = abs2(T̃ - ω * c / a)
    B = (a * b - abs2(c)) / a^2
    Δ² = A + B * (ω^2 + a * bath.W^2)
    Δ² * Gₕ(ω, bath, a)
end


"""
    function SA_marcus(t_idx, tf, TG, C, bath::HybridOhmicBath, i = 1, j =2)

Calculate the Marcus rate of Smirnov-Amin ME -- ``Γ_{ij}``.
"""
function SA_marcus(t_idx, tf, TG, C, bath::HybridOhmicBath, i = 1, j = 2)
    ω = 2π * (TG.ω[t_idx, j] - TG.ω[t_idx, i])
    a = C.a[t_idx, j, i]
    b = C.b[t_idx, j, i]
    c = C.c[t_idx, j, i]
    d = C.d[t_idx, j, i]
    T = TG.T[t_idx, j, i]
    G = TG.G[t_idx, j, i]

    T̃ = T - 1.0im * G / tf - d * bath.ϵl
    A = abs2(T̃ - ω * c / a)
    B = (b - abs2(c) / a) * bath.W^2
    Δ² = A + B
    Δ² * Gₗ(ω, bath, a)
end


"""
    function SA_Γ(t_idx, tf, TG, C, bath::HybridOhmicBath, i = 1, j =2)

Calculate the hybrid rate of Smirnov-Amin ME -- ``Γ_{ij}``.
"""
function SA_Γ(t_idx, tf, TG, C, bath::HybridOhmicBath, i = 1, j = 2)
    ω = 2π * (TG.ω[t_idx, j] - TG.ω[t_idx, i])
    a = C.a[t_idx, j, i]
    b = C.b[t_idx, j, i]
    c = C.c[t_idx, j, i]
    d = C.d[t_idx, j, i]
    T = TG.T[t_idx, j, i]
    G = TG.G[t_idx, j, i]

    T̃ = T - 1.0im * G / tf - d * bath.ϵ
    A = abs2(T̃ - ω * c / a)
    B = (a * b - abs2(c)) / a^2
    integrand(x) = (A + B * (x^2 + a * bath.W^2)) * Gₗ(ω - x, bath, a) * Gₕ(x, bath, a)

    # split the integration range for better accuracy
    μ = ω - a * bath.ϵl
    σ = sqrt(a) * bath.width_l
    γ = a * bath.width_h

    integration_range = sort([Inf, -Inf, μ - 3*σ, μ, μ + 3*σ, -3*γ, 0, 3*γ, 2*σ, -2*σ])

    res, err = quadgk(integrand, integration_range..., rtol = 1e-7, atol = 1e-7)
    res/2/π, err/2/π
end


function SA_τ(t_idx, tf, TG, C, bath::HybridOhmicBath, i = 1, j = 2)
    ω = 2π * (TG.ω[t_idx, j] - TG.ω[t_idx, i])
    a = C.a[t_idx, j, i]
    1 / max(abs(ω), sqrt(a)*bath.W)
end


function SA_ΓMatrix(tf, TG, C, bath::HybridOhmicBath)
    t_dim, lvl = size(TG.ω)
    Γ_data = Array{Float64, 3}(undef, lvl-1, lvl, t_dim)
    for j = 1:lvl
        for i = 1:(j-1)
            for t_idx = 1:t_dim
                Γ_data[i, j, t_idx] = SA_Γ(t_idx, tf, TG, C, bath, i, j)[1]
            end
        end
        for i = j:(lvl-1)
            for t_idx = 1:t_dim
                Γ_data[i, j, t_idx] = SA_Γ(t_idx, tf, TG, C, bath, i+1, j)[1]
            end
        end
    end
    ΓMatrix(range(0,1,length=t_dim), Γ_data)
end


"""
    function solve_SA(TG, C, bath::HybridOhmicBath, u0, tf; kwargs...)

Solve Smirnov-Amin ME. Currently this function is experimental. It does not have the same interface as the other solvers.
"""
function solve_SA(TG, C, bath::HybridOhmicBath, u0, tf; kwargs...)
    Γ = SA_ΓMatrix(tf, TG, C, bath::HybridOhmicBath)
    prob = ODEProblem{true}(Γ, u0, (0.0, 1.0), tf)
    solve(prob; alg_hints = [:nonstiff], kwargs...)
end


function SA_lz_rotate(sol, θ_itp, lvl; s = nothing)
    if s == nothing
        s = sol.t
    end
    y = Matrix{Float64}(undef, length(s), length(lvl))
    for (i, v) in enumerate(s)
        U = QTBase.@unitary_landau_zener(θ_itp(v))
        ρ = U*diagm(sol(v))*U'
        y[i, :] = diag(ρ)[lvl]
    end
    y
end
