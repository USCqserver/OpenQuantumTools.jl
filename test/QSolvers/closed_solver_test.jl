using OpenQuantumTools, Test, OrdinaryDiffEq

H = DenseHamiltonian([(s) -> 1.0], -[σz], unit = :ħ)
u0 = PauliVec[1][1]
annealing = Annealing(H, u0)

tf = π
U = exp(1.0im * tf * σz)
sol = solve_schrodinger(annealing, tf, alg = Tsit5(), reltol = 1e-4)
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

sol = solve_von_neumann(annealing, tf, alg = Tsit5(), reltol = 1e-4)
@test sol(tf) ≈ U * u0 * u0' * U' atol = 1e-4 rtol = 1e-4

sol = solve_von_neumann(
    annealing,
    tf,
    alg = Exprb32(),
    vectorize = true,
    reltol = 1e-4,
)
@test reshape(sol(tf), 2, 2) ≈ U * u0 * u0' * U' atol = 1e-4 rtol = 1e-4

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
sol = solve_unitary(annealing, tf, alg = Tsit5(), callback = cb, retol = 1e-4)
@test sol(tf) ≈ σx atol = 1e-4 rtol = 1e-4
cb = InstPulseCallback([0.5 * tf], (c, i) -> c .= σx * c * σx)
sol =
    solve_von_neumann(annealing, tf, alg = Tsit5(), callback = cb, retol = 1e-4)
@test sol(tf) ≈ u0 * u0' rtol = 1e-4 atol = 1e-4

# test for SparseHamiltonian
Hp = spσz ⊗ spσz
H = SparseHamiltonian([(s) -> 1.0], -[Hp], unit = :ħ)
u0 = PauliVec[1][1] ⊗ PauliVec[1][1]
annealing = Annealing(H, u0)

tf = 1.0
U = exp(1.0im * tf * Array(Hp))
sol = solve_schrodinger(annealing, tf, alg = Tsit5(), retol=1e-4)
@test sol(1.0) ≈ U * u0 atol = 1e-4 rtol = 1e-4
sol = solve_schrodinger(annealing, tf, alg = Exprb32(), retol=1e-4)
@test sol(1.0) ≈ U * u0 atol = 1e-4 rtol = 1e-4
@test_broken solve_schrodinger(
    annealing,
    tf,
    alg = TRBDF2(),
    abstol = 1e-8,
    reltol = 1e-8,
)
#@test_broken sol(1.0) ≈ U * u0 atol = 1e-4 rtol = 1e-4 # TRBDF2() cannot correctly start
sol = solve_unitary(annealing, tf, alg = Tsit5(), retol=1e-4)
@test sol(tf) ≈ U atol = 1e-4 rtol = 1e-4
sol = solve_von_neumann(annealing, tf, alg = Tsit5(), retol=1e-4)
@test sol(tf) ≈ U * u0 * u0' * U' atol = 1e-4 rtol = 1e-4
sol = solve_von_neumann(annealing, tf, alg = MagnusGauss4(), vectorize=true, dt=1/50)
@test sol[end] ≈ (U * u0 * u0' * U')[:] atol = 1e-4 rtol = 1e-4
