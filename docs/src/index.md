# HOME

*This is a quantum toolbox for Julia programming language.*

This package has the following features

## Math Symbols
The package defines some commonly used mathematical symbols and operators. Some of its features include:
### Pauli Matrices
Single qubit Pauli matrices are defined by ```σx```, ```σy```, ```σz``` and ```σi```. One can simply calculate the tensor product of single qubits Pauli matrices by
```julia-repl
julia> σx⊗σz
4×4 Array{Complex{Float64},2}:
 0.0+0.0im   0.0+0.0im  1.0+0.0im   0.0+0.0im
 0.0+0.0im  -0.0+0.0im  0.0+0.0im  -1.0+0.0im
 1.0+0.0im   0.0+0.0im  0.0+0.0im   0.0+0.0im
 0.0+0.0im  -1.0+0.0im  0.0+0.0im  -0.0+0.0im
```
The eigenvectors of each Pauli matrices are also defined in constant [`PauliVec`](@ref), where ```PauliVec[1]```,```PauliVec[2]```,```PauliVec[3]``` corresponds to eigenvectors of ```σx```, ```σy```, ```σz```. Additionally, the first eigenvector is the one with positive eigenvalue
```julia-repl
julia> σz*PauliVec[3][1] == PauliVec[3][1]
true
```
Additionally, a sparse version of Pauli matrices are defined in ```spσx``` et al. They can be used to construct Hamiltonian in the form of sparse matrix.
### Construction of Multi-Qubits Matrices
Multi-qubits matrices can be construct with two different methods.

## Hamiltonian
