module QHamil

import LinearAlgebra:kron, tr, mul!, axpy!, ishermitian, Hermitian, eigmin, I, eigen
import LinearAlgebra.BLAS:her!
import SparseArrays:issparse
import Arpack:eigs
import Optim:optimize # optimize package is used to find minimal gap

include("../QTBase/QTBase.jl")
import .QTBase:σx, σz, σy, σi, σ, ⊗, PauliVec, comm, comm!, spσx, spσz, spσi, spσy

export matrix_decompose, check_positivity

export q_translate, construct_hamming_weight_op, ising_terms, standard_driver, collective_operator, GHZ_entanglement_witness, local_field_term, two_local_term

export Hamiltonian, inst_population, gibbs_state, eigen_eval, eigen_state_continuation!, low_level_hamiltonian, minimum_gap

include("matrix_tool.jl")
include("utility_func.jl")
include("hamil_obj.jl")

end  # module QHamil
