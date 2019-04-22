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
    res = zeros(eltype(h.m[1]), size(h.m[1]))
    for (f,m) in zip(h.f,h.m)
        axpy!(f(t),m,res)
    end
    res
end

function Hamiltonian(f::Array{Function,1},m::Array{Array{T,2},1}) where T<:Number
    Hamiltonian(f,m,log2(size(m[1], 1)))
end

struct UnitlessAdiabaticFrameHamiltonian
    dθfun::Array{Function, 1}
    hfun::Array{Function, 1}
    dθmat::Array{Array{T, 2}, 1} where T<: Number
    hmat::Array{Array{T, 2}, 1} where T<: Number
    n_qubit::Int64
    _tf_cache::Array{Float64, 1}
end

struct UnitAdiabaticFrameHamiltonian
    dθfun::Array{Function, 1}
    hfun::Array{Function, 1}
    dθmat::Array{Array{T, 2}, 1} where T<: Number
    hmat::Array{Array{T, 2}, 1} where T<: Number
    n_qubit::Int64
    _tf_cache::Array{Float64, 1}
end

function AdiabaticFrameHamiltonian(dθfun, hfun, dθmat, hmat; unitless=true)
    if unitless == true
        return UnitlessAdiabaticFrameHamiltonian(dθfun, hfun, dθmat, hmat, log2(size(hmat[1], 1)), [1.0])
    else
        return UnitAdiabaticFrameHamiltonian(dθfun, hfun, dθmat, hmat, log2(size(hmat[1], 1)), [1.0])
    end
end

function set_tf!(H::UnitlessAdiabaticFrameHamiltonian, tf)
    scale = tf / H._tf_cache[1]
    H._tf_cache[1] = tf
    for i in eachindex(H.hmat)
        H.hmat[i] *= scale
    end
end

function set_tf!(H::UnitAdiabaticFrameHamiltonian, tf)
    scale = tf / H._tf_cache[1]
    H._tf_cache[1] = tf
    for i in eachindex(H.dθmat)
        H.dθmat[i] /= scale
    end
end

function (h::UnitlessAdiabaticFrameHamiltonian)(s::Real)
    res = zeros(eltype(h.hmat[1]), size(h.hmat[1]))
    for (f,m) in zip(h.dθfun, h.dθmat)
        axpy!(f(s), m, res)
    end
    for (f,m) in zip(h.hfun, h.hmat)
        axpy!(f(s), m, res)
    end
    res
end

function (h::UnitAdiabaticFrameHamiltonian)(t::Real)
    s = t / h._tf_cache[1]
    res = zeros(eltype(h.hmat[1]), size(h.hmat[1]))
    for (f,m) in zip(h.dθfun, h.dθmat)
        axpy!(f(s), m, res)
    end
    for (f,m) in zip(h.hfun, h.hmat)
        axpy!(f(s), m, res)
    end
    res
end

struct UnitlessAdiabaticFramePausingHamiltonian
    dθfun::Array{Function, 1}
    hfun::Array{Function, 1}
    dθmat::Array{Array{T, 2}, 1} where T<: Number
    hmat::Array{Array{T, 2}, 1} where T<: Number
    sp::Float64
    sd::Float64
    tf_ext::Float64
    n_qubit::Int64
end

struct UnitAdiabaticFramePausingHamiltonian
    dθfun::Array{Function, 1}
    hfun::Array{Function, 1}
    dθmat::Array{Array{T, 2}, 1} where T<: Number
    hmat::Array{Array{T, 2}, 1} where T<: Number
    sp::Float64
    sd::Float64
    tf::Float64
    tf_ext::Float64
    n_qubit::Int64
end


function construct_pausing_hamiltonian(sp, sd, H::UnitlessAdiabaticFrameHamiltonian)
    tf_ext = (1+sd)*H._tf_cache[1]
    hmat = [(1+sd)*x for x in H.hmat]
    UnitlessAdiabaticFramePausingHamiltonian(H.dθfun, H.hfun, H.dθmat, hmat, sp, sd, tf_ext, H.n_qubit)
end

function construct_pausing_hamiltonian(sp, sd, H::UnitAdiabaticFrameHamiltonian)
    tf_ext = (1+sd)*H._tf_cache[1]
    dθmat = [x/(1+sd) for x in H.dθmat]
    UnitAdiabaticFramePausingHamiltonian(H.dθfun, H.hfun, dθmat, H.hmat, sp, sd, H._tf_cache[1], tf_ext, H.n_qubit)
end


function (h::UnitlessAdiabaticFramePausingHamiltonian)(s::Real)
    res = zeros(eltype(h.hmat[1]), size(h.hmat[1]))
    if s <= h.sp
        for (f,m) in zip(h.dθfun, h.dθmat)
            axpy!(f(s), m, res)
        end
        for (f,m) in zip(h.hfun, h.hmat)
            axpy!(f(s), m, res)
        end
    elseif s <= h.sp + h.sd
        for (f,m) in zip(h.hfun, h.hmat)
            axpy!(f(h.sp), m, res)
        end
    else
        sn = s-h.sd
        for (f,m) in zip(h.dθfun, h.dθmat)
            axpy!(f(sn), m, res)
        end
        for (f,m) in zip(h.hfun, h.hmat)
            axpy!(f(sn), m, res)
        end
    end
    res
end

function (h::UnitAdiabaticFramePausingHamiltonian)(t::Real)
    s = t / h.tf
    res = zeros(eltype(h.hmat[1]), size(h.hmat[1]))
    if s <= h.sp
        for (f,m) in zip(h.dθfun, h.dθmat)
            axpy!(f(s), m, res)
        end
        for (f,m) in zip(h.hfun, h.hmat)
            axpy!(f(s), m, res)
        end
    elseif s <= h.sp + h.sd
        for (f,m) in zip(h.hfun, h.hmat)
            axpy!(f(h.sp), m, res)
        end
    else
        sn = s-h.sd
        for (f,m) in zip(h.dθfun, h.dθmat)
            axpy!(f(sn), m, res)
        end
        for (f,m) in zip(h.hfun, h.hmat)
            axpy!(f(sn), m, res)
        end
    end
    res
end
