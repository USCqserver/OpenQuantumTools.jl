using QuantumAnnealingTools, Test, Random
using OrdinaryDiffEq

H = DenseHamiltonian([(s) -> 1.0], -[σz], unit=:ħ)
u0 = PauliVec[1][1]
coupling = ConstantCouplings(["Z"], unit=:ħ)
bath = EnsembleFluctuator([1.0, 2.0], [2.0, 1.0])
annealing = Annealing(H, u0, coupling=coupling, bath=bath)
tf = 2.0

prob = build_ensembles(annealing, tf, :stochastic)
Random.seed!(1234)
sol1 = solve(prob, Tsit5(), EnsembleSerial(), trajectories=1, save_everystep=false)

Random.seed!(1234)
sol2 = solve(prob, Tsit5(), EnsembleSerial(), trajectories=1, save_everystep=false)

@test sol1[1][end] ≈ sol2[1][end] atol = 1e-6 rtol = 1e-6

sol3 = solve(prob, Tsit5(), EnsembleSerial(), trajectories=1, save_everystep=false)

@test !(isapprox(sol1[1][end], sol3[1][end], atol=1e-6, rtol=1e-6))