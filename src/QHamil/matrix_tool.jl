const σx = [0.0+0.0im 1; 1 0]
const σy = [0.0+0.0im -1.0im; 1.0im 0]
const σz = [1.0+0.0im 0; 0 -1]
const σi = [1.0+0.0im 0; 0 1]
const σ = [σx, σy, σz, σi]
xvec = [[1.0+0.0im, 1.0]/sqrt(2), [1.0+0.0im, -1.0]/sqrt(2)]
yvec = [[1.0im, -1.0]/sqrt(2), [1.0im, 1.0]/sqrt(2)]
zvec = [[1.0+0.0im, 0], [0, 1.0+0.0im]]

const spσz = sparse(σz)
const spσx = sparse(σx)
const spσy = sparse(σy)
const spσi = sparse(I, 2, 2)

"""
    PauliVec

Constants for the eigenvectors of single qubit Pauli matrices. Indices 1, 2 and 3 corresponds to the eigenvectors of ``σ_x``, ``σ_y`` and ``σ_z``.

# Examples
```julia-repl
julia> σx*PauliVec[1][1] == PauliVec[1][1]
true
```
"""
const PauliVec = [xvec, yvec, zvec]

"""
    ⊗(A, B)

Calculate the tensor product of `A` and `B`.

# Examples
```julia-repl
julia> σx⊗σz
4×4 Array{Complex{Float64},2}:
 0.0+0.0im   0.0+0.0im  1.0+0.0im   0.0+0.0im
 0.0+0.0im  -0.0+0.0im  0.0+0.0im  -1.0+0.0im
 1.0+0.0im   0.0+0.0im  0.0+0.0im   0.0+0.0im
 0.0+0.0im  -1.0+0.0im  0.0+0.0im  -0.0+0.0im
```
"""
⊗ = kron

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
function q_translate(h::String; sp=false)
    # define operator map to replace [XYZI] with corresponding Pauli matrices
    function operator_map(x)
        if sp ==false
            σ_tag = "σ"
        else
            σ_tag = "spσ"
        end
        res = ""
        for i in range(1, length=length(x)-1)
            res = res * σ_tag * lowercase(x[i]) * "⊗"
        end
        res = res * σ_tag * lowercase(x[end])
    end

    h_str = replace(h, r"[XYZI]+" => operator_map)
    eval(Meta.parse(h_str))
end

"""
    ising_terms(ops, q_ind, weight, num_qubit; sp=false)

Construct an ising term of the multi-qubits Hamiltonian. `ops` is a list of Pauli operator names which appears in the ising term. `q_ind` is the list of indices corresponding to the Pauli matrices in `ops`. `weight` is the constant factor of this ising term. `num_qubit` is the total number of qubits. A sparse matrix can be construct by setting `sp` to `true`. The following example construct an ising term of `` Z_1 I Z_3/2 ``.

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
function ising_terms(ops, q_ind, weight, num_qubit; sp=false)
    if sp == false
        σ_tag = "σ"
        i_tag = σi
    else
        σ_tag = "spσ"
        i_tag = spσi
    end
    res = weight
    for i in 1:num_qubit
        idx = findfirst((x)->x==i, q_ind)
        if idx != nothing
            op2 = eval(Meta.parse(σ_tag * lowercase(ops[idx])))
            res = res ⊗ op2
        else
            res = res ⊗ i_tag
        end
    end
    res
end

"""
    collective_operator(op, num_qubit; sp=false)

Construct the collective operator for a system of `num_qubit` qubits. `op` is the name of the collective Pauli matrix. For example, the following code construct an `` IZ + ZI `` matrix. Generate sparse matrix when `sp` is set to true.

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
function collective_operator(op, num_qubit; sp=false)
    op_name = uppercase(op)
    res = ""
    for idx in 1:num_qubit
        res = res * "I"^(idx-1)*op_name*"I"^(num_qubit-idx) *"+"
    end
    q_translate(res[1:end-1]; sp=sp)
end

"""
    standard_driver(num_qubit; sp=false)

Construct the standard driver Hamiltonian for a system of `num_qubit` qubits. For example, a two qubits standard driver matrix is `` IX + XI ``. Generate sparse matrix when sp is set to true.
"""
function standard_driver(num_qubit; sp=false)
    res = ""
    for idx in 1:num_qubit
        res = res * "I"^(idx-1)*"X"*"I"^(num_qubit-idx) *"+"
    end
    q_translate(res[1:end-1], sp=sp)
end

"""
    construct_hamming_weight_op(num_qubit::Int64, op::String; sp=false)

Construct the Hamming weight operator for system of size `num_qubit`. The type of the Hamming weight operator is specified by op: "x", "y" or "z". Generate sparse matrix when sp is set to true.

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
function construct_hamming_weight_op(num_qubit::Int64, op::String; sp=false)
    0.5 * (num_qubit*I - collective_operator(op, num_qubit=num_qubit, sp=sp))
end

function GHZ_entanglement_witness(num_qubit)
    s = collective_operator("z", num_qubit)
    for k in 2:num_qubit
        s += ising_terms(["z","z"], [k-1,k], 1.0, num_qubit)
    end
    (num_qubit-1)I - s
end

function local_field_term(h, idx, num_qubit; sp=false)
    res = ising_terms(["z"], [idx[1]], h[1], num_qubit, sp=sp)
    for i in 2:length(idx)
        res += ising_terms(["z"], [idx[i]], h[i], num_qubit, sp=sp)
    end
    res
end

function two_local_term(j, idx, num_qubit; sp=false)
    res = ising_terms(["z", "z"], idx[1], j[1], num_qubit, sp=sp)
    for i in 2:length(idx)
        res += ising_terms(["z", "z"], idx[i], j[i], num_qubit, sp=sp)
    end
    res
end
