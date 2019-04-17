function solve_adiabatic_me(hfun, u0, tf, inter_op, bath)

end

function adiabatic_me_update!(du, u, A, γ, S)
    A2 = abs2.(A)
    γA = γ .* A2
    Γ = sum(γA, dims=1)
    dim = size(du)[1]
    for a in 1:dim
        for b in 1:a-1
            du[a, a] += γA[a, b] * u[b, b] - γA[b, a] * u[a, a]
            du[a, b] += -0.5 * (Γ[a] + Γ[b]) * u[a, b] + γ[1, 1] * A[a, a] * A[b, b] * u[a, b]
        end
        for b in a+1:dim
            du[a, a] += γA[a, b] * u[b, b] - γA[b, a] * u[a, a]
            du[a, b] += -0.5 * (Γ[a] + Γ[b]) * u[a, b] + γ[1, 1] * A[a, a] * A[b, b] * u[a, b]
        end
    end
    H_ls = Diagonal(sum(S .* A2, dims=1)[1,:])
    axpy!(-1.0im, H_ls*u-u*H_ls, du)
end
