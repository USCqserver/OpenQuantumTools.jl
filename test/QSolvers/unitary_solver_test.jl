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

tf = 2.0
sol = solve_unitary(annealing, tf, alg = Tsit5())
@test sol(1.0) ≈ exp(1.0im * tf * σz) atol = 1e-4 rtol = 1e-4
sol = solve_unitary(annealing, tf, alg = Tsit5(), dimensionless_time = false)
@test sol(tf) ≈ exp(1.0im * tf * σz) atol = 1e-4 rtol = 1e-4
sol =
    solve_unitary(annealing, tf, alg = TRBDF2(), vectorize = true, abstol = 1e-10, reltol = 1e-10)
@test reshape(sol(1.0), 2, 2) ≈ exp(1.0im * tf * σz) atol = 1e-4 rtol = 1e-4

control = InstPulseControl([0.5], (x) -> σx)
annealing = Annealing(H, u0, control = control)
sol = solve_unitary(annealing, tf, alg = Tsit5())
@test sol(1.0) ≈ exp(0.5im * tf * σz) * σx * exp(0.5im * tf * σz) atol = 1e-5
sol = solve_unitary(annealing, tf, alg = Tsit5(), dimensionless_time = false)
@test sol(2.0) ≈ exp(0.5im * tf * σz) * σx * exp(0.5im * tf * σz) atol = 1e-5

control = InstDEPulseControl([0.5], (x) -> σx, :state)
annealing = Annealing(H, u0, control = control)
sol = solve_unitary(
    annealing,
    tf,
    alg = Tsit5(),
    de_array_constructor = de_wrapper,
)
@test sol(1.0) ≈ exp(0.5im * tf * σz) * σx * exp(0.5im * tf * σz) atol = 1e-5
@test sol.u[end].state == 2
