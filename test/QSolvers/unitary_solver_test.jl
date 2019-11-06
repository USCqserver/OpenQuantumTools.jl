using QuantumAnnealingTools, Test, OrdinaryDiffEq

H = DenseHamiltonian([(s)->1.0], -[σz], unit=:ħ)
u0 = PauliVec[1][1]
annealing = Annealing(H, u0)

tf = 1.0
sol = solve_unitary(annealing, tf, alg=Tsit5())
@test sol(1.0) ≈ exp(1.0im*tf*σz) atol = 1e-4 rtol = 1e-4
