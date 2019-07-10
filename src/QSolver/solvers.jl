function load_diff_eq()
    @eval QuantumAnnealingTools begin
        using DifferentialEquations:ODEProblem, solve, Tsit5, ODEFunction
    end
end

function solve_schrodinger(𝐇, u0; kwargs...)
    function f(du, u, p, t)
        hmat = -1.0im * 𝐇(t)
        mul!(du, hmat, u)
    end
    function f_jac(J, u, p, t)
        hmat = 𝐇(t)
        mul!(J, -1.0im, hmat)
    end
    ff =ODEFunction(f; jac=f_jac)
    prob = ODEProblem(ff, u0, (0.0, 1.0))
    sol = solve(prob, Tsit5(); kwargs...)
end

function solve_von_neumann(𝐇, u0; kwargs...)
    function f(du, u, p ,t)
        hmat = -1.0im * 𝐇(t)
        mul!(du, hmat, u)
        axpy!(-1.0, u*hmat, du)
    end
    prob = ODEProblem(f, u0, (0.0,1.0))
    sol = solve(prob,Tsit5(); kwargs...)
end

function adiabatic_frame_ame(hfun, u0, inter_op, γf, sf; rtol=1e-6, atol=1e-6)
    function f(du, u, p, t)
        hmat = -1.0im * hfun(t)
        mul!(du, hmat, u)
        axpy!(-1.0, u*hmat, du)
        ω = diag(hmat)
        ω_ba = repeat(ω, 1, length(ω))
        ω_ba = transpose(ω_ba) - ω_ba
        γm = p*γf.(ω_ba)
        sm = p*sf.(ω_ba)
        for op in inter_op
            adiabatic_me_update!(du, u, op(t), γm, sm)
        end
    end
    prob = ODEProblem(f, u0, (0.0, 1.0), tf)
    sol = solve(prob, Tsit5(), reltol=rtol, abstol=atol, save_everystep=false)
end

function solve_adiabatic_me(hfun, u0, tf, inter_op, γf, sf; rtol=1e-6, atol=1e-6)
    function f(du, u, p ,t)
        hmat = hfun(t)
        w,v = eigen!(Hermitian(hmat))
        ρ = v' * u * v
        H = Diagonal(w)
        dρ = -1.0im * p * (H * ρ - ρ * H)
        ω_ba = repeat(w, 1, length(w))
        ω_ba = transpose(ω_ba) - ω_ba
        γm = p*γf.(ω_ba)
        sm = p*sf.(ω_ba)
        A = v' * inter_op * v
        adiabatic_me_update!(dρ, ρ, A, γm, sm)
        mul!(du, v, dρ*v')
    end
    prob = ODEProblem(f, u0, (0.0,1.0), tf)
    sol = solve(prob, Tsit5(), reltol=rtol, abstol=atol)
end
