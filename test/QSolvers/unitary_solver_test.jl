using QuantumAnnealingTools, Test, OrdinaryDiffEq

H = DenseHamiltonian([(s)->1.0], -[σz], unit=:ħ)
u0 = PauliVec[1][1]
annealing = Annealing(H, u0)

tf = 2.0
sol = solve_unitary(annealing, tf, alg=Tsit5())
@test sol(1.0) ≈ exp(1.0im*tf*σz) atol = 1e-4 rtol = 1e-4
sol = solve_unitary(annealing, tf, alg=TRBDF2(), abstol=1e-10, reltol=1e-10)
@test sol(1.0) ≈ exp(1.0im*tf*σz) atol = 1e-4 rtol = 1e-4

control = InstPulseControl([0.5], (x)->σx)
annealing = Annealing(H, u0, control=control)
sol = solve_unitary(annealing, tf, alg=Tsit5())
@test sol(1.0) ≈ exp(0.5im*tf*σz)*σx*exp(0.5im*tf*σz) atol=1e-5
sol = solve_unitary(annealing, tf, alg=Tsit5(), span_unit=true)
@test sol(2.0) ≈ exp(0.5im*tf*σz)*σx*exp(0.5im*tf*σz) atol=1e-5
