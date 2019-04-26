module QTBase

using Reexport
import LinearAlgebra:kron, mul!, axpy!, I, ishermitian, Hermitian, eigmin, eigen, tr, eigen!
import LinearAlgebra.BLAS:her!
import SparseArrays:sparse, issparse, spzeros, SparseMatrixCSC
import Arpack:eigs
import Optim:optimize

export temperature_2_beta, temperature_2_freq, beta_2_temperature, freq_2_temperature

export σx, σz, σy, σi, σ, ⊗, PauliVec, spσx, spσz, spσi, spσy

export q_translate, construct_hamming_weight_op, ising_terms, standard_driver, collective_operator, GHZ_entanglement_witness, local_field_term, two_local_term

export matrix_decompose, check_positivity

export inst_population, gibbs_state, eigen_eval, eigen_state_continuation!, low_level_hamiltonian, minimum_gap, proj_low_lvl

include("unit_util.jl")
include("math_util.jl")
include("matrix_util.jl")
include("proj_util.jl")


include("../QInterpolate/QInterpolate.jl")
include("../Integration/Integration.jl")
@reexport using .QInterpolate
@reexport using .Integration

end  # module QTBase
