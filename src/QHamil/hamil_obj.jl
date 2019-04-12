"""
Object to hold general time dependent Hamiltonians, whose form is assumed to be a summation of time dependent function times constant matrices.

**Fields**
- `f` -- list of time dependent functions.
- `m` -- list of constant matrices.
- `n_qubit` -- total number of qubits.
"""
struct Hamiltonian
    " List of time dependent functions "
    f::Array{Function, 1}
    " List of constant matrices "
    m::Array{Array{ComplexF64,2},1}
    " Total number of qubits "
    n_qubit::Int64
end

" Evaluate the Hamiltonian at time `t` "
function (h::Hamiltonian)(t::Real)
    res = zeros(ComplexF64, size(h.m[1]))
    for (f,m) in zip(h.f,h.m)
        axpy!(f(t),m,res)
    end
    res
end

function Hamiltonian(f::Array{Function,1},m::Array{Array{T,2},1}) where T<:Number
    Hamiltonian(f,m,log2(size(m[1])[1]))
end
