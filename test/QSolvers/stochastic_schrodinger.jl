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
sol1 = solve(prob, Tsit5(), EnsembleSerial(), trajectories=1)
tstops1 = []
for i = 2:length(sol1[1].t)
    if sol1[1].t[i] == sol1[1].t[i - 1]
        push!(tstops1, sol1[1].t[i])
    end
end

Random.seed!(1234)
sol2 = solve(prob, Tsit5(), EnsembleSerial(), trajectories=1)
tstops2 = []
for i = 2:length(sol2[1].t)
    if sol2[1].t[i] == sol2[1].t[i - 1]
        push!(tstops2, sol2[1].t[i])
    end
end

@test sol1[1][end] ≈ sol2[1][end] atol = 1e-6 rtol = 1e-6
@test !isempty(tstops1)
@test tstops1 == tstops2
