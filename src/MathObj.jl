export σx, σz, σy, σi, ⊗, plus, minus, comm, comm!, Hamiltonian

const σx = [0.0+0.0im 1; 1 0]
const σy = [0.0+0.0im -1.0im; 1.0im 0]
const σz = [1.0+0.0im 0; 0 -1]
const σi = [1.0+0.0im 0; 0 1]
const plus = [1.0+0.0im, 1.0] / sqrt(2)
const minus = [1.0+0.0im, -1.0] / sqrt(2)

⊗ = kron

"""
    comm(A, B)

Calculate the commutator of matrices 'A' and 'B'
"""
function comm(A, B)
    A*B - B*A
end

"""
    comm!(Y, A, B)

Calculate the commutator of matrices 'A' and 'B'. Write the result in 'Y'
"""
function comm!(Y, A, B)
    mul!(Y,A,B)
    axpy!(-1.0,B*A,Y)
end

" Object to hold general time dependent Hamiltonians, whose form is assumed to be a summation of time dependent function times constant matrices. "
struct Hamiltonian
    " List of time dependent functions "
    f::Array{Function, 1}
    " List of constant matrices "
    m::Array{Array{ComplexF64,2},1}
    " Total number of qubits "
    n_qubit::Int64
end

" Evaluate the Hamiltonian at time 't' "
function (h::Hamiltonian)(t::Real)
    res = zeros(h.m[1])
    for (f,m) in zip(h.f,h.m)
        axpy!(f(t),m,res)
    end
    res
end

function Hamiltonian(f::Array{Function,1},m::Array{Array{T,2},1}) where T<:Number
    Hamiltonian(f,m,log2(size(m[1])[1]))
end
