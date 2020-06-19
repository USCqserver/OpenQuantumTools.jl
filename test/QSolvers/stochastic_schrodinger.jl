using QuantumAnnealingTools, Test, Random
using OrdinaryDiffEq

mutable struct DENoiseArray{T,N} <: DEDataArray{T,N}
    x::Array{T,N}
    n::Vector{Float64}
end
de_wrapper(x) = DENoiseArray(x, [1.0])

control = InstPulseControl([0.5], (x) -> σx)

H = DenseHamiltonian([(s) -> 1.0], -[σz], unit = :ħ)
u0 = PauliVec[1][1]
coupling = ConstantCouplings(["Z"], unit = :ħ)
bath = EnsembleFluctuator([1.0, 2.0], [2.0, 1.0])
annealing = Annealing(H, u0, coupling = coupling, bath = bath)
tf = 2.0

Random.seed!(1234)
exp_res = [
    0.6026564365923035 + 0.373882986790114im,
    0.6026564365923035 - 0.373882986790114im,
]
exp_jump_t = [
    0.17221764281694701
    0.276846651563355
    0.2822012637167232
    0.30558258606841276
    0.37969915405033106
    0.4592929326126119
    0.4934539602929163
    0.6349091730846184
]

sol = solve_stochastic_schrodinger(annealing, tf, alg = Tsit5())
@test sol[end] == exp_res

Random.seed!(1234)
sol = solve_stochastic_schrodinger(
    annealing,
    tf,
    alg = Tsit5(),
    de_array_constructor = de_wrapper,
    fluctuator_de_field = :n,
)
@test sol[end] == exp_res
left_n = [s.n for s in sol(exp_jump_t, continuity=:left)]
right_n = [s.n for s in sol(exp_jump_t, continuity=:right)]
@test !any((x)->x[1]==x[2], zip(left_n, right_n))

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
