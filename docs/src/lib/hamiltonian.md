# Hamiltonian
## Hamiltonian types
```@docs
DenseHamiltonian
DenseHamiltonian(funcs, mats; unit = :h, EIGS = EIGEN_DEFAULT)
SparseHamiltonian
SparseHamiltonian(funcs, mats; unit = :h, EIGS = EIGEN_DEFAULT)
AdiabaticFrameHamiltonian
AdiabaticFrameHamiltonian(ωfuns, geofuns)
```
##
```@docs
evaluate(H::OpenQuantumBase.AbstractHamiltonian, s::Real)
eigen_decomp(H::OpenQuantumBase.AbstractHamiltonian, s; lvl::Int=2)
eigen_decomp(H::OpenQuantumBase.AbstractHamiltonian, s::AbstractArray{Float64,1}; lvl::Int=2)
```