# Matrix Utilities
## Symbols
```@docs
PauliVec
⊗
```
## Math utilities
```@docs
check_positivity(m::AbstractMatrix)
check_unitary(𝐔::AbstractMatrix; rtol = 1e-6, atol = 1e-8)
check_density_matrix(ρ; atol::Real=0, rtol::Real=atol>0 ? 0 : √eps())
fidelity(ρ, σ)
partial_trace(ρ::Matrix, qubit_2_keep)
matrix_decompose(mat::AbstractMatrix, basis::Vector{<:AbstractMatrix})
```
## Construction utilities
```@docs
q_translate(h::String; sp = false)
q_translate_state(h::String; normal=false)
single_clause(ops, q_ind, weight, num_qubit; sp=false)
standard_driver(num_qubit; sp = false)
local_field_term(h, idx, num_qubit; sp=false)
two_local_term(j, idx, num_qubit; sp=false)
collective_operator(op, num_qubit; sp=false)
```