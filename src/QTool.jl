module QTool

using LinearAlgebra
include("MathObj.jl")
include("Unitary.jl")
include("QUnit.jl")
export q_translate, construct_hamming_weight_op

# == translate string to hamiltonian matrices ==
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
    construct_hamming_weight_op(num_qubit::Int64, op::String)

Construct the Hamming weight operator for system of size 'num_qubit'. The type of the Hamming weight operator is specified by op: "x", "y" or "z".

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

end # end module
