using QuantumAnnealingTools, Test, OrdinaryDiffEq

# test InstDEPulseControl
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

sol = solve_von_neumann(annealing, tf, alg = Tsit5())
@test sol(1.0) ≈ U * u0 * u0' * U' atol = 1e-4 rtol = 1e-4
sol = solve_von_neumann(
    annealing,
    tf,
    alg = TRBDF2(),
    abstol = 1e-10,
    reltol = 1e-10,
)
@test sol(1.0) ≈ U * u0 * u0' * U' atol = 1e-4 rtol = 1e-4
sol = solve_von_neumann(annealing, tf, alg = Exprb32(), vectorize = true)
@test reshape(sol(1.0), 2, 2) ≈ U * u0 * u0' * U' atol = 1e-4 rtol = 1e-4
sol =
    solve_von_neumann(annealing, tf, alg = Tsit5(), dimensionless_time = false)
@test sol(tf) ≈ U * u0 * u0' * U' atol = 1e-4 rtol = 1e-4

control = InstPulseControl([0.5], (x) -> σx)
annealing = Annealing(H, u0, control = control)
sol = solve_von_neumann(annealing, tf, alg = Tsit5())
@test sol[end] ≈ u0 * u0' atol = 1e-4 rtol = 1e-4
sol = solve_von_neumann(annealing, tf, alg = Tsit5(), vectorize = true)
@test sol[end] ≈ (u0*u0')[:] atol = 1e-4 rtol = 1e-4


control = InstDEPulseControl([0.5], (x) -> σx, :state)
annealing = Annealing(H, u0, control = control)
sol = solve_von_neumann(
    annealing,
    tf,
    alg = Tsit5(),
    de_array_constructor = de_wrapper,
)
@test sol.u[end] ≈ u0 * u0' atol = 1e-4 rtol = 1e-4
@test sol.u[end].state == 2

sol = solve_von_neumann(
    annealing,
    tf,
    alg = Tsit5(),
    de_array_constructor = de_wrapper,
    vectorize = true,
)
@test sol[end] ≈ (u0*u0')[:] atol = 1e-4 rtol = 1e-4

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
