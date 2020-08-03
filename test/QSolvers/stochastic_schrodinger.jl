using QuantumAnnealingTools, Test, Random
using OrdinaryDiffEq

H = DenseHamiltonian([(s) -> 1.0], -[σz], unit = :ħ)
u0 = PauliVec[1][1]
coupling = ConstantCouplings(["Z"], unit = :ħ)
bath = EnsembleFluctuator([1.0, 2.0], [2.0, 1.0])
annealing = Annealing(H, u0, coupling = coupling, bath = bath)
tf = 2.0

prob = build_ensembles(annealing,tf,:stochastic)
Random.seed!(1234)
sol = solve(prob, Tsit5(), EnsembleSerial(), trajectories=1)
tstops = []
for i = 2:length(sol[1].t)
    if sol[1].t[i] == sol[1].t[i-1]
        push!(tstops, sol[1].t[i])
    end
end
exp_res = [
    0.6595793462134859 + 0.2607910764094174im
    0.6595793462134859 - 0.2607910764094174im
]
exp_jump_t = [
    0.30221824355127086
    0.9554780593261047
    1.2999133449599989
    1.5091713624528147
    1.5198805867595513
    1.5666432314629304
    1.7148763674267669
    1.8740639245513286
    1.9423859799119376
]

@test sol[1][end] == exp_res
@test tstops == exp_jump_t
