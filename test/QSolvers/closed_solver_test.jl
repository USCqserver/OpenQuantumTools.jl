using OpenQuantumTools, Test, OrdinaryDiffEq

H = DenseHamiltonian([(s) -> 1.0], -[σz], unit = :ħ)
u0 = PauliVec[1][1]
ρ₀ = 0.2 * PauliVec[1][1] * PauliVec[1][1]' + 0.8 * PauliVec[1][2] * PauliVec[1][2]'
annealing = Annealing(H, u0)
annealingρ = Annealing(H, ρ₀)

tf = π
U = exp(1.0im * tf * σz)
sol = solve_schrodinger(annealing, tf, alg = Tsit5(), reltol = 1e-4)
@test_throws DiffEqBase.NoDefaultAlgorithmError solve_schrodinger(annealing, tf, reltol = 1e-4)
@test sol(tf) ≈ U * u0 atol = 1e-4 rtol = 1e-4
sol = solve_schrodinger(annealing, tf, alg = Exprb32(), reltol = 1e-4)
@test sol(tf) ≈ U * u0 atol = 1e-4 rtol = 1e-4
# TRBDF2 did not perform well for this case
sol = solve_schrodinger(
    annealing,
    tf,
    alg = TRBDF2(),
    abstol = 1e-9,
    reltol = 1e-9,
)
@test sol(tf) ≈ U * u0 atol = 1e-4 rtol = 1e-4

sol = solve_unitary(annealing, tf, alg = Tsit5(), reltol = 1e-4)
@test sol(tf) ≈ U atol = 1e-4 rtol = 1e-4
sol = solve_unitary(
    annealing,
    tf,
    alg = TRBDF2(),
    reltol = 1e-9,
    abstol = 1e-9,
    vectorize = true,
)
@test sol(tf) ≈ U[:] atol = 1e-4 rtol = 1e-4

@test_logs (:warn, "The initial state is a pure state. It is more efficient to use the Schrodinger equation solver.") solve_von_neumann(annealing, tf, alg = Tsit5(), reltol = 1e-4)
sol = solve_von_neumann(annealing, tf, alg = Tsit5(), reltol = 1e-4)
@test sol(tf) ≈ U * u0 * u0' * U' atol = 1e-4 rtol = 1e-4

sol = solve_von_neumann(
    annealingρ,
    tf,
    alg = Exprb32(),
    vectorize = true,
    reltol = 1e-4,
)
@test reshape(sol(tf), 2, 2) ≈ U * ρ₀' * U' atol = 1e-4 rtol = 1e-4

# test control protocols on Schrodinger solver
# test InstPulseControl
cb = InstPulseCallback([0.5 * tf], (c, i) -> c .= σx * c)
sol = solve_schrodinger(
    annealing,
    tf,
    alg = Tsit5(),
    callback = cb,
    reltol = 1e-4,
)
@test sol.u[end] ≈ u0 atol = 1e-4 rtol = 1e-4
sol = solve_unitary(annealing, tf, alg = Tsit5(), callback = cb, reltol = 1e-4)
@test sol(tf) ≈ σx atol = 1e-4 rtol = 1e-4
cb = InstPulseCallback([0.5 * tf], (c, i) -> c .= σx * c * σx)
sol =
    solve_von_neumann(annealingρ, tf, alg = Tsit5(), callback = cb, reltol = 1e-4)
@test sol(tf) ≈ ρ₀ rtol = 1e-4 atol = 1e-4

# test for SparseHamiltonian
Hp = spσz ⊗ spσz
H = SparseHamiltonian([(s) -> 1.0], -[Hp], unit = :ħ)
u0 = PauliVec[1][1] ⊗ PauliVec[1][1]
u1 = PauliVec[1][2] ⊗ PauliVec[1][2]
ρ₀ = 0.2*u0*u0' + 0.8*u1*u1' 
annealing = Annealing(H, u0)
annealingρ = Annealing(H, ρ₀)

tf = 1.0
U = exp(1.0im * tf * Array(Hp))
sol = solve_schrodinger(annealing, tf, alg = Tsit5(), reltol=1e-4)
@test sol(1.0) ≈ U * u0 atol = 1e-4 rtol = 1e-4
sol = solve_schrodinger(annealing, tf, alg = Exprb32(), reltol=1e-4)
@test sol(1.0) ≈ U * u0 atol = 1e-4 rtol = 1e-4
@test_broken solve_schrodinger(
    annealing,
    tf,
    alg = TRBDF2(),
    abstol = 1e-8,
    reltol = 1e-8,
)
#@test_broken sol(1.0) ≈ U * u0 atol = 1e-4 rtol = 1e-4 # TRBDF2() cannot correctly start
sol = solve_unitary(annealing, tf, alg = Tsit5(), reltol=1e-4)
@test sol(tf) ≈ U atol = 1e-4 rtol = 1e-4
sol = solve_von_neumann(annealingρ, tf, alg = Tsit5(), reltol=1e-4)
@test sol(tf) ≈ U * ρ₀ * U' atol = 1e-4 rtol = 1e-4
sol = solve_von_neumann(annealingρ, tf, alg = MagnusGauss4(), vectorize=true, dt=1/50)
@test sol[end] ≈ (U * ρ₀ * U')[:] atol = 1e-4 rtol = 1e-4

f(s) = -σz
H = hamiltonian_from_function(f)
u0 = PauliVec[1][1]
ρ₀ = 0.2 * PauliVec[1][1] * PauliVec[1][1]' + 0.8 * PauliVec[1][2] * PauliVec[1][2]'
annealing = Annealing(H, u0)
annealingρ = Annealing(H, ρ₀)
tf = π
sol = solve_schrodinger(annealing, tf, alg = Tsit5(), reltol = 1e-4)
U = exp(1.0im * tf * σz)
@test sol(tf) ≈ U * u0 atol = 1e-4 rtol = 1e-4
sol = solve_von_neumann(annealing, tf, alg = Tsit5(), reltol = 1e-4)
@test sol(tf) ≈ U * u0 * u0' * U' atol = 1e-4 rtol = 1e-4
sol = solve_unitary(annealing, tf, alg = Tsit5(), reltol = 1e-4)
@test sol(tf) ≈ U atol = 1e-4 rtol = 1e-4