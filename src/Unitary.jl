using DifferentialEquations
export calculate_unitary

function calculate_unitary(𝐇;rtol=1e-6,atol=1e-8)
    u0 = Matrix{ComplexF64}(I, size(𝐇(0.0)))
    function f(du, u, p, t)
        hmat = -1.0im * 𝐇(t)
        mul!(du, hmat, u)
    end
    prob = ODEProblem(f, u0, (0.0,1.0))
    sol = solve(prob,Tsit5(),reltol=rtol,abstol=atol)
end
