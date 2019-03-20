module QTool

using Reexport
using Optim # optimize package is used to find minimal gap
@reexport using LinearAlgebra
@reexport using DifferentialEquations
include("hamiltonian_construction_util.jl")
include("hamiltonian_obj.jl")
include("unit_util.jl")
include("unitary_util.jl")
include("math_util.jl")



export q_translate, construct_hamming_weight_op, ising_terms, standard_driver, collective_operator

export ħ, Planck, Boltzmann
export temperature_2_beta, temperature_2_freq

export calculate_unitary, unitary_check, solve_schrodinger, solve_von_neumann

export Hamiltonian, eigen_value_eval, eigen_state_eval, inst_population, gibbs_state, eigen_sys_eval, eigen_state_continuation!, low_level_hamiltonian, minimum_gap


end # end module
