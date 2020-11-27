# Math Symbols
`OpenQuantumTools` defines some commonly used mathematical symbols and operators, including:

## Pauli matrices
The constants `σx`, `σy`, `σz` and `σi` represent the Pauli matrices. The binary operator `⊗` (can be typed as `\otimes<tab>`) represents the tensor product. For example,
```julia-repl
julia> σx⊗σz
4×4 Array{Complex{Float64},2}:
 0.0+0.0im   0.0+0.0im  1.0+0.0im   0.0+0.0im
 0.0+0.0im  -0.0+0.0im  0.0+0.0im  -1.0+0.0im
 1.0+0.0im   0.0+0.0im  0.0+0.0im   0.0+0.0im
 0.0+0.0im  -1.0+0.0im  0.0+0.0im  -0.0+0.0im
```
calculates the tensor product of `σx` and `σz`.

The eigenvectors of each Pauli matrix are also stored in the constant [`PauliVec`](@ref), where `PauliVec[1]`,`PauliVec[2]`,`PauliVec[3]` correspond to the  eigenvectors of `σx`, `σy`, `σz` respectively. The first element in each `PauliVec[i]` is the one with a positive eigenvalue
```julia-repl
julia> σz*PauliVec[3][1] == PauliVec[3][1]
true
```
Additionally, sparse versions of the Pauli matrices are defined in ```spσx```, ```spσy```, ```spσz``` and ```spσi```. They can be used to construct `SparseHamiltonian`.

## Utility functions
`OpenQuantumTools` provides various utility functions for convenience.

* [`check_positivity`](@ref): check if a matrix is positive semi-definite.
* [`check_unitary`](@ref): check if a matrix is unitary.
* [`check_density_matrix`](@ref): check if a matrix is a valid density matrix.
* [`fidelity`](@ref): calculate the fidelity between two density matrices `ρ` and `σ` using ``Tr[\sqrt{\sqrt{\rho}\sigma\sqrt{\rho}}]^2``

```julia-repl
julia> ρ = PauliVec[1][1]*PauliVec[1][1]'
julia> σ = PauliVec[3][1]*PauliVec[3][1]'
julia> fidelity(ρ, σ)
0.49999999999999944
```

* [`partial_trace`](@ref): calculate the partial trace of a matrix

```julia-repl
julia> ρ1 = [0.4 0.2; 0.2 0.6]; ρ2 = [0.5 0; 0 0.5];
julia> partial_trace(ρ1⊗ρ2, [1])
2×2 Array{Float64,2}:
0.4  0.2
0.2  0.6
```

* [`matrix_decompose`](@ref): project a matrix onto a list of basis elements

```julia-repl
julia> matrix_decompose(1.0*σx+2.0*σy+3.0*σz, [σx,σy,σz])
3-element Array{Complex{Float64},1}:
1.0 + 0.0im
2.0 + 0.0im
3.0 + 0.0im
```

## Construction of multi-qubit matrices
A collection of functions to conveniently construct multi-qubit matrices is also listed below, to which a keyword argument `sp` can be supplied to generate sparse matrices.

* [`standard_driver`](@ref) builds the standard driver Hamiltonian in quantum annealing:
```julia-repl
julia> standard_driver(2) == σx⊗σi + σi⊗σx
true
julia> standard_driver(2, sp=true) == spσx⊗spσi + spσi⊗spσx
true
```

*  [`q_translate`](@ref) builds a multi-qubit matrix from its string representation:
```julia-repl
julia> q_translate("ZZ+0.5ZI-XZ")
4×4 Array{Complex{Float64},2}:
  1.5+0.0im   0.0+0.0im  -1.0+0.0im   0.0+0.0im
  0.0+0.0im  -0.5+0.0im   0.0+0.0im   1.0+0.0im
 -1.0+0.0im   0.0+0.0im  -1.5+0.0im  -0.0+0.0im
  0.0+0.0im   1.0+0.0im  -0.0+0.0im   0.5+0.0im
```

* [`collective_operator`](@ref) constructs a collective Pauli operator (the same Pauli operator acting on each individual qubit):
```julia-repl
julia> collective_operator("z", 3) == σz⊗σi⊗σi + σi⊗σz⊗σi + σi⊗σi⊗σz
true
```

* [`single_clause`] builds a single-clause term:
```julia-repl
julia> single_clause(["z","z"], [2,3], -2, 4) == -2σi⊗σz⊗σz⊗σi
true
```

* [`local_field_term`](@ref) builds a local field terms of the form ``ΣhᵢZᵢ``:
```julia-repl
julia> local_field_term([1.0, 0.5], [1, 2], 2) == σz⊗σi+0.5σi⊗σz
true
```

* [`two_local_term`](@ref) builds two-local terms of the form ``∑JᵢⱼZᵢZⱼ``:
```julia-repl
julia> two_local_term([1.0, 0.5], [[1,2], [1,3]], 3) == σz⊗σz⊗σi + 0.5σz⊗σi⊗σz
true
```

## Construction of multi-qubit states
The quantum state of a spin system can be constructed by[`q_translate_state`](@ref)
```julia-repl
julia> q_translate_state("001")
8-element Array{Complex{Float64},1}:
 0.0 + 0.0im
 1.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im
 0.0 + 0.0im
```
In the string representation, `0` and `1` represent the eigenstates of the ``σ_z`` operator.
