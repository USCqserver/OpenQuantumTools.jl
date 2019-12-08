using QuantumAnnealingTools, Test, OrdinaryDiffEq

H = DenseHamiltonian([(s) -> 1.0], -[σz], unit = :ħ)
u0 = PauliVec[1][1]
annealing = Annealing(H, u0)

tf = 1.0
sol = solve_von_neumann(annealing, tf, alg = Tsit5())
U = exp(1.0im * tf * σz)
@test sol(1.0) ≈ U * u0 * u0' * U' atol = 1e-4 rtol = 1e-4
sol = solve_von_neumann(annealing, tf, alg = TRBDF2(), abstol = 1e-8, reltol = 1e-8)
@test sol(1.0) ≈ U * u0 * u0' * U' atol = 1e-4 rtol = 1e-4
sol = solve_von_neumann(annealing, tf, alg = Exprb32(), vectorize = true)
@test reshape(sol(1.0), 2, 2) ≈ U * u0 * u0' * U' atol = 1e-4 rtol = 1e-4


# test for sparse matrix
Hp = spσz ⊗ spσz
H = SparseHamiltonian([(s) -> 1.0], -[Hp], unit = :ħ)
u0 = PauliVec[1][1] ⊗ PauliVec[1][1]
annealing = Annealing(H, u0)

tf = 1.0
U = exp(1.0im * tf * Array(Hp))

sol = solve_von_neumann(annealing, tf, alg = Exprb32(), vectorize = true)
@test reshape(sol(1.0), 4, 4) ≈ U * u0 * u0' * U'

sol = solve_von_neumann(annealing, tf, alg = Tsit5())
@test sol(1.0) ≈ U * u0 * u0' * U' atol = 1e-4 rtol = 1e-4
