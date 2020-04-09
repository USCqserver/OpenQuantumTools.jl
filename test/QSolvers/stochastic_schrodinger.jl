using QuantumAnnealingTools, Test
using OrdinaryDiffEq

mutable struct DENoiseArray{T,N} <: DEDataArray{T,N}
    x::Array{T,N}
    n::Vector{Float64}
end
de_wrapper(x) = DENoiseArray(x, [1.0])

control = InstPulseControl([0.5], (x) -> σx)

H = DenseHamiltonian([(s)->1.0], -[σz], unit=:ħ)
u0 = PauliVec[1][1]
coupling = ConstantCouplings(["Z"], unit=:ħ)
bath = EnsembleFluctuator([1.0, 2.0], [2.0, 1.0])
annealing = Annealing(H, u0, coupling=coupling, bath=bath)
tf = 2.0

#sol = solve_stochastic_schrodinger(annealing, tf, alg=Tsit5())
#sol = solve_stochastic_schrodinger(annealing, tf, alg=Tsit5(), de_array_constructor=de_wrapper, fluctuator_de_field=:n)

# prob, callback = build_ensemble_problem(annealing, tf, :stochastic_schrodinger)
# sol = solve(prob, Tsit5(), EnsembleSerial(); trajectories=2, callback=callback, abstol=1e-6, reltol=1e-6);

#prob, callback = build_ensemble_problem(annealing, tf, :stochastic_schrodinger, de_array_constructor=de_wrapper, fluctuator_de_field=:n)

# annealing = Annealing(H, u0, coupling=coupling, bath=bath, control=control)
# sol = solve_stochastic_schrodinger(annealing, 2.0, alg=Tsit5())


# mutable struct DEHybridArray{T,N} <: DEDataArray{T,N}
#     x::Array{T,N}
#     n::Vector{Float64}
#     s::Int
# end
# de_wrapper(x) = DEHybridArray(x, [1.0], 1)
#
# control = InstDEPulseControl([0.5], (x) -> σx, :s)
# annealing = Annealing(H, u0, coupling=coupling, bath=bath, control=control)
#
# sol = solve_stochastic_schrodinger(annealing, 2.0, alg=Tsit5(), de_array_constructor=de_wrapper, fluctuator_de_field=:n)

#prob, callback = build_ensemble_problem(annealing, tf, :stochastic_schrodinger, de_array_constructor=de_wrapper, additional_field=:n)
