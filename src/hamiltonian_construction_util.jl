"""
    q_translate(h::String)

Convert a string `h` representing multi-qubits Pauli matrices summation into its numerical form.

# Examples
```julia-repl
julia> q_translate("X+2.0Z")
2×2 Array{Complex{Float64},2}:
 2.0+0.0im   1.0+0.0im
 1.0+0.0im  -2.0+0.0im
```
"""
function q_translate(h::String)
    h_str = replace(h, r"[XYZI]+" => operator_map)
    @debug "Translated string" h_str
    eval(Meta.parse(h_str))
end

function operator_map(x)
    res = ""
    for i in range(1, length=length(x)-1)
        res = res * "σ"*lowercase(x[i]) * "⊗"
    end
    res = res * "σ"*lowercase(x[end])
end

"""
    ising_terms(ops, q_ind, weight, num_qubit)

Construct an ising term of the multi-qubits Hamiltonian. `ops` is a list of Pauli operator names which appears in the ising term. `q_ind` is the list of indices corresponding to the Pauli matrices in `ops`. `weight` is the constant factor of this ising term. `num_qubit` is the total number of qubits. For example, an ising term `` Z_1 I Z_3/2 `` can be constructed as

# Examples
```julia-repl
julia> ising_terms(["z", "z"], [1, 3], 0.5, 3)
8×8 Array{Complex{Float64},2}:
 0.5+0.0im   0.0+0.0im  0.0+0.0im   0.0+0.0im   0.0+0.0im   0.0+0.0im   0.0+0.0im   0.0+0.0im
 0.0+0.0im  -0.5+0.0im  0.0+0.0im  -0.0+0.0im   0.0+0.0im  -0.0+0.0im   0.0+0.0im  -0.0+0.0im
 0.0+0.0im   0.0+0.0im  0.5+0.0im   0.0+0.0im   0.0+0.0im   0.0+0.0im   0.0+0.0im   0.0+0.0im
 0.0+0.0im  -0.0+0.0im  0.0+0.0im  -0.5+0.0im   0.0+0.0im  -0.0+0.0im   0.0+0.0im  -0.0+0.0im
 0.0+0.0im   0.0+0.0im  0.0+0.0im   0.0+0.0im  -0.5+0.0im  -0.0+0.0im  -0.0+0.0im  -0.0+0.0im
 0.0+0.0im  -0.0+0.0im  0.0+0.0im  -0.0+0.0im  -0.0+0.0im   0.5-0.0im  -0.0+0.0im   0.0-0.0im
 0.0+0.0im   0.0+0.0im  0.0+0.0im   0.0+0.0im  -0.0+0.0im  -0.0+0.0im  -0.5+0.0im  -0.0+0.0im
 0.0+0.0im  -0.0+0.0im  0.0+0.0im  -0.0+0.0im  -0.0+0.0im   0.0-0.0im  -0.0+0.0im   0.5-0.0im
```
"""
function ising_terms(ops, q_ind, weight, num_qubit)
    res = weight
    for i in 1:num_qubit
        idx = findfirst((x)->x==i, q_ind)
        if idx != nothing
            op2 = eval(Meta.parse("σ"*lowercase(ops[idx])))
            res = res ⊗ op2
        else
            res = res ⊗ σi
        end
    end
    res
end

"""
    collective_operator(op, num_qubit)

Construct the collective operator for a system of `num_qubit` qubits. `op` is the name of the collective Pauli matrix. For example, the following code construct an `` IZ + ZI `` matrix.

# Examples
```julia-repl
julia> collective_operator("z", 2)
4×4 Array{Complex{Float64},2}:
 2.0+0.0im  0.0+0.0im  0.0+0.0im   0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im   0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im   0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  -2.0+0.0im
```
"""
function collective_operator(op, num_qubit)
    op_name = uppercase(op)
    res = ""
    for idx in 1:num_qubit
        res = res * "I"^(idx-1)*op_name*"I"^(num_qubit-idx) *"+"
    end
    q_translate(res[1:end-1])
end

"""
    standard_driver(num_qubit)

Construct the standard driver Hamiltonian for a system of `num_qubit` qubits. For example, a two qubits standard driver matrix is `` IX + XI ``.
"""
function standard_driver(num_qubit)
    res = ""
    for idx in 1:num_qubit
        res = res * "I"^(idx-1)*"X"*"I"^(num_qubit-idx) *"+"
    end
    q_translate(res[1:end-1])
end

"""
    construct_hamming_weight_op(num_qubit::Int64, op::String)

Construct the Hamming weight operator for system of size `num_qubit`. The type of the Hamming weight operator is specified by op: "x", "y" or "z".

# Examples
```julia-repl
julia> construct_hamming_weight_op(2,"z")
4×4 Array{Complex{Float64},2}:
 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  1.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  1.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  2.0+0.0im
```
"""
function construct_hamming_weight_op(num_qubit::Int64, op::String)
    res = zeros(ComplexF64, (2^num_qubit, 2^num_qubit))
    for i in range(1, length=num_qubit)
        op_name = "I"^(i-1) * op * "I"^(num_qubit-i)
        operator = I - eval(Meta.parse(operator_map(op_name)))
        res = res + operator
    end
    rmul!(res, 0.5)
end
