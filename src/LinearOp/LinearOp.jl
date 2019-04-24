module LinearOp

import LinearAlgebra:axpy!
import LinearAlgebra.BLAS.gemm!

export LinearOperator, update!, comm!, AdiabaticFrameHamiltonian
export set_tf!, construct_pausing_hamiltonian
export UnitlessAdiabaticFrameHamiltonian, UnitAdiabaticFrameHamiltonian, UnitlessAdiabaticFramePausingHamiltonian, UnitAdiabaticFramePausingHamiltonian

include("linear_operator.jl")
include("hamil_obj.jl")

end  # module QHamil
