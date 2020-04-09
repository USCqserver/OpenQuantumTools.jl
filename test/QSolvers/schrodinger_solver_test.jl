using QuantumAnnealingTools, Test, OrdinaryDiffEq

# test structures for InstDEPulseControl
mutable struct DEStateMachineArray{T,N} <: DEDataArray{T,N}
    """Array data"""
    x::Array{T,N}
    """Current state"""
    state::Int
end
de_wrapper(x) = DEStateMachineArray(x, 1)

H = DenseHamiltonian([(s) -> 1.0], -[σz], unit = :ħ)
u0 = PauliVec[1][1]
annealing = Annealing(H, u0)

tf = π
U = exp(1.0im * tf * σz)
sol = solve_schrodinger(annealing, tf, alg = Tsit5())
@test sol(1.0) ≈ U * u0 atol = 1e-4 rtol = 1e-4
sol = solve_schrodinger(annealing, tf, alg = Exprb32())
@test sol(1.0) ≈ U * u0 atol = 1e-4 rtol = 1e-4
# TRBDF2 did not perform well for this case
sol = solve_schrodinger(
    annealing,
    tf,
    alg = TRBDF2(),
    abstol = 1e-9,
    reltol = 1e-9,
)
@test sol(1.0) ≈ U * u0 atol = 1e-4 rtol = 1e-4
sol =
    solve_schrodinger(annealing, tf, alg = Tsit5(), dimensionless_time = false)
@test sol.t[end] ≈ tf
@test sol.u[end] ≈ U * u0 atol = 1e-4 rtol = 1e-4

# test control protocols on Schrodinger solver
# test InstPulseControl
control = InstPulseControl([0.5], (x) -> σx)
annealing = Annealing(H, u0, control = control)
sol = solve_schrodinger(annealing, tf, alg = Tsit5())
@test sol.u[end] ≈ u0 atol = 1e-4 rtol = 1e-4
sol =
    solve_schrodinger(annealing, tf, alg = Tsit5(), dimensionless_time = false)
@test sol.u[end] ≈ u0 atol = 1e-4 rtol = 1e-4

# test InstDEPulseControl
control = InstDEPulseControl([0.5], (x) -> σx, :state)
annealing = Annealing(H, u0, control = control)
sol = solve_schrodinger(
    annealing,
    tf,
    alg = Tsit5(),
    de_array_constructor = de_wrapper,
)

@test sol.u[end] ≈ u0 atol = 1e-4 rtol = 1e-4
@test sol.u[end].state == 2

@test_throws ErrorException solve_schrodinger(annealing, tf, alg = Tsit5())

sol = solve_schrodinger(
    annealing,
    tf,
    alg = Tsit5(),
    de_array_constructor = de_wrapper,
    dimensionless_time = false,
)

@test sol(tf) ≈ u0 atol = 1e-4 rtol = 1e-4
@test sol.u[end].state == 2


# test for SparseHamiltonian
Hp = spσz ⊗ spσz
H = SparseHamiltonian([(s) -> 1.0], -[Hp], unit = :ħ)
u0 = PauliVec[1][1] ⊗ PauliVec[1][1]
annealing = Annealing(H, u0)

tf = 1.0
sol = solve_schrodinger(annealing, tf, alg = Tsit5())
@test sol(1.0) ≈ exp(1.0im * tf * Array(Hp)) * u0 atol = 1e-4 rtol = 1e-4
sol = solve_schrodinger(annealing, tf, alg = Exprb32())
@test sol(1.0) ≈ exp(1.0im * tf * Array(Hp)) * u0 atol = 1e-4 rtol = 1e-4
sol = solve_schrodinger(
    annealing,
    tf,
    alg = TRBDF2(),
    abstol = 1e-8,
    reltol = 1e-8,
)
@test sol(1.0) ≈ exp(1.0im * tf * Array(Hp)) * u0 atol = 1e-4 rtol = 1e-4
