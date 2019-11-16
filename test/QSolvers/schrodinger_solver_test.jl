using QuantumAnnealingTools, Test, OrdinaryDiffEq

H = DenseHamiltonian([(s)->1.0], -[σz], unit=:ħ)
u0 = PauliVec[1][1]
annealing = Annealing(H, u0)

tf = 1.0
sol = solve_schrodinger(annealing, tf, alg=Tsit5())
@test sol(1.0) ≈ exp(1.0im*tf*σz)*u0 atol = 1e-4 rtol = 1e-4
sol = solve_schrodinger(annealing, tf, alg=Exprb32())
@test sol(1.0) ≈ exp(1.0im*tf*σz)*u0 atol = 1e-4 rtol = 1e-4
# TRBDF2 did not perform well for this case
sol = solve_schrodinger(annealing, tf, alg=TRBDF2(), abstol=1e-8, reltol=1e-8)
@test sol(1.0) ≈ exp(1.0im*tf*σz)*u0 atol = 1e-4 rtol = 1e-4

tf = π
sol = solve_schrodinger(annealing, tf, alg=Tsit5(), span_unit=true)
@test sol.t[end] ≈ tf
@test sol.u[end] ≈ exp(1.0im*tf*σz)*u0 atol = 1e-4 rtol = 1e-4

# test the parallel solver
sol = solve_schrodinger(annealing, [1.0, π], Tsit5(), span_unit=true, tstops=[0.5])
@test 0.5 in sol[1].t
@test sol[1](1.0) ≈ exp(1.0im*σz)*u0 atol = 1e-4 rtol = 1e-4
@test 0.5*π in sol[2].t
@test sol[2].u[end] ≈ exp(1.0im*π*σz)*u0 atol = 1e-4 rtol = 1e-4

sol = solve_schrodinger(annealing, [1.0, π], Tsit5(), EnsembleThreads(), span_unit=true, tstops=[0.5])
@test 0.5 in sol[1].t
@test sol[2].t[end] ≈ π

Hp = spσz⊗spσz
H = SparseHamiltonian([(s)->1.0], -[Hp], unit=:ħ)
u0 = PauliVec[1][1]⊗PauliVec[1][1]
annealing = Annealing(H, u0)

tf = 1.0
sol = solve_schrodinger(annealing, tf, alg=Tsit5())
@test sol(1.0) ≈ exp(1.0im*tf*Array(Hp))*u0 atol = 1e-4 rtol = 1e-4
sol = solve_schrodinger(annealing, tf, alg=Exprb32())
@test sol(1.0) ≈ exp(1.0im*tf*Array(Hp))*u0 atol = 1e-4 rtol = 1e-4
sol = solve_schrodinger(annealing, tf, alg=TRBDF2(), abstol=1e-8, reltol=1e-8)
@test sol(1.0) ≈ exp(1.0im*tf*Array(Hp))*u0 atol = 1e-4 rtol = 1e-4
