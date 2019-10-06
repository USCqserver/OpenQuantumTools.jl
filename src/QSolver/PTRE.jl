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
