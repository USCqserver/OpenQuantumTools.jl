# Matrix Utilities
## Symbols
```@docs
PauliVec
⊗
```
## Construction Utilities
```@docs
q_translate(h::String; sp = false)
q_translate_state(h::String; normal=false)
single_clause(ops, q_ind, weight, num_qubit; sp=false)
standard_driver(num_qubit; sp = false)
local_field_term(h, idx, num_qubit; sp=false)
two_local_term(j, idx, num_qubit; sp=false)
collective_operator(op, num_qubit; sp=false)
```