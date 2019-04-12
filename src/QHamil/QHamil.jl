module QHamil

import LinearAlgebra:kron, tr, mul!, axpy!, ishermitian, Hermitian, eigmin, I, eigen
import LinearAlgebra.BLAS:her!
import SparseArrays:sparse
import Arpack:eigs
import Optim:optimize # optimize package is used to find minimal gap

export σx, σz, σy, σi, σ, ⊗, PauliVec, comm, comm!, spσx, spσz, spσi, spσy

export matrix_decompose, check_positivity

export q_translate, construct_hamming_weight_op, ising_terms, standard_driver, collective_operator, GHZ_entanglement_witness, local_field_term, two_local_term

export Hamiltonian, eigen_value_eval, eigen_state_eval, inst_population, gibbs_state, eigen_sys_eval, eigen_state_continuation!, low_level_hamiltonian, minimum_gap, sp_eigen_sys_eval

include("matrix_tool.jl")
include("utility_func.jl")
include("hamil_obj.jl")



"""
    comm(A, B)

Calculate the commutator of matrices `A` and `B`.
"""
function comm(A, B)
    A*B - B*A
end

"""
    comm!(Y, A, B)

Calculate the commutator of matrices `A` and `B` and add the result to `Y`.
"""
function comm!(Y, A, B)
    mul!(Y,A,B)
    axpy!(-1.0,B*A,Y)
end


end  # module QHamil
