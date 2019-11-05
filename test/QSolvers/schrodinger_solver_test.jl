using QuantumAnnealingTools, Test, OrdinaryDiffEq

H = DenseHamiltonian([(s)->1.0], -[σz], unit=:ħ)
u0 = PauliVec[1][1]
annealing = Annealing(H, u0)

tf = π
sol = solve_schrodinger(annealing, tf, alg=Tsit5())
@test sol(1.0) ≈ exp(-1.0im*tf*σz)*u0 atol = 1e-4 rtol = 1e-4

sol = solve_schrodinger(annealing, tf, alg=Tsit5(), span_unit=true)
@test sol.t[end] ≈ tf
@test sol.u[end] ≈ exp(-1.0im*tf*σz)*u0 atol = 1e-4 rtol = 1e-4
