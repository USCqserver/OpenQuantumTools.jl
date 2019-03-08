# HOME

*This is a quantum toolbox for Julia programming language.*

This package has the following features

## Define mathematical symbols and operations
The package defines some commonly used mathematical symbols and operators. For example, the single qubit Pauli matrices are defined by ```σx```, ```σy```, ```σz``` and ```σi```. One can simply calculate the tensor product of single qubits Pauli matrices by
```julia-repl
julia> σx⊗σz
4×4 Array{Complex{Float64},2}:
 0.0+0.0im   0.0+0.0im  1.0+0.0im   0.0+0.0im
 0.0+0.0im  -0.0+0.0im  0.0+0.0im  -1.0+0.0im
 1.0+0.0im   0.0+0.0im  0.0+0.0im   0.0+0.0im
 0.0+0.0im  -1.0+0.0im  0.0+0.0im  -0.0+0.0im
```
