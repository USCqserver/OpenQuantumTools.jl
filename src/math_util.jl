"""
    minimum_gap(h)

Calculate the minimum gap of Hamiltonian `h` using Optim.jl package.
"""
function minimum_gap(h)
    function gap(s)
        eig_sys = eigen(Hermitian(h(s)))
        eig_sys.values[2]-eig_sys.values[1]
    end
    optimize(gap, 0.0, 1.0)
end


pulse_on_density!(c::Vector, p) = c .= conj(p) ⊗ p * c
pulse_on_density!(c::Matrix, p) = c .= p * c * p'
