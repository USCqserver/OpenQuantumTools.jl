using QuantumAnnealingTools, Test

load_diff_eq()
hfun(t) = 5*σx
sol = calculate_unitary(hfun)
u_res = exp(-1.0im*5*0.5*σx)
@test isapprox(u_res, sol(0.5), rtol=1e-6, atol=1e-8)
@test check_unitary(u_res)
@test !check_unitary([0 1; 0 0])
ode_res = solve_schrodinger(hfun, PauliVec[3][1])
@test isapprox(exp(-1.0im*5*σx) * PauliVec[3][1] , ode_res.u[end], rtol=1e-4, atol=1e-4)
