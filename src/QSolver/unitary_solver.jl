function solve_unitary(A::Annealing, tf::Real; kwargs...)
    u0 = Matrix{ComplexF64}(I, size(A.H(0.0)))
    p = AnnealingParams(A.H, float(tf))
    ff = ODEFunction(mul_ode; jac = mul_jac)
    prob = ODEProblem(ff, u0, A.sspan, p)
    # currently stiff algorithm does not support complex type
    solve(prob; alg_hints = [:nonstiff], tstops=A.tstops, kwargs...)
end

function uni_f(du, u, p, t)
    p.H(tf, t)
end
