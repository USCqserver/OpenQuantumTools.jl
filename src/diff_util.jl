function adiabatic_me_update!(du, u, A, γ, S)
    dim = size(du)[1]
    for a in 1:dim
        for b in 1:dim
            A = abs2(A[a,b])
            du[a,a] += γ[a,b]*A*u[b,b]
            du[b,b] -= γ[a,b]*A*u[b,b]
            if a!=b
                du[a,b] -= (0.5*A*(γ[a,b]+γ[b,a]) + 0.5*γ[1,1]*(A[a,a]+A[b,b])^2)*u[a,b]
            end
        end
    end
    H_ls = Diagonal(sum(S .* abs2(A), dims=1))
    axpy!(-1.0im, H_ls*u-u*H_ls, du)
end
