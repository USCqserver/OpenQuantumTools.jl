module QTool

using Reexport
using Optim # optimize package is used to find minimal gap
using Arpack
@reexport using SparseArrays
@reexport using LinearAlgebra

# functions needed for noise module
import SpecialFunctions.trigamma

# end import for noise module

include("load_util.jl")
include("hamiltonian_construction_util.jl")
include("hamiltonian_obj.jl")
include("unit_util.jl")
include("unitary_util.jl")
include("math_util.jl")
include("noise_util.jl")

export load_diff_eq

export σx, σz, σy, σi, σ, ⊗, PauliVec, comm, comm!, spσx, spσz, spσi, spσy
export matrix_decompose, check_positivity

export q_translate, construct_hamming_weight_op, ising_terms, standard_driver, collective_operator, GHZ_entanglement_witness, local_field_term, two_local_term

export ħ, Planck, Boltzmann
export temperature_2_beta, temperature_2_freq

export calculate_unitary, unitary_check, solve_schrodinger, solve_von_neumann

export Hamiltonian, eigen_value_eval, eigen_state_eval, inst_population, gibbs_state, eigen_sys_eval, eigen_state_continuation!, low_level_hamiltonian, minimum_gap, sp_eigen_sys_eval

export OhmicBath, γ, correlation, HybridOhmic


end # end module
