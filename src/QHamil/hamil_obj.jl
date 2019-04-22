struct UnitlessAdiabaticFrameHamiltonian
    dθ::LinearOperator
    h::LinearOperator
    n_qubit::Int64
    _tf_cache::Array{Float64, 1}
end

struct UnitAdiabaticFrameHamiltonian
    dθ::LinearOperator
    h::LinearOperator
    n_qubit::Int64
    _tf_cache::Array{Float64, 1}
end

function AdiabaticFrameHamiltonian(dθfun, hfun, dθmat, hmat; unitless=true)
    dθ_op = LinearOperator(dθfun, dθmat)
    h_op = LinearOperator(hfun, hmat)
    n_qubit = log2(size(hmat[1], 1))
    if unitless == true
        return UnitlessAdiabaticFrameHamiltonian(dθ_op, h_op, n_qubit, [1.0])
    else
        return UnitAdiabaticFrameHamiltonian(dθ_op, h_op, n_qubit, [1.0])
    end
end

function set_tf!(H::UnitlessAdiabaticFrameHamiltonian, tf)
    scale = tf / H._tf_cache[1]
    H._tf_cache[1] = tf
    multiply!(H.h, scale)
end

function set_tf!(H::UnitAdiabaticFrameHamiltonian, tf)
    scale = tf / H._tf_cache[1]
    H._tf_cache[1] = tf
    multiply!(H.dθ, 1/scale)
end

function (h::UnitlessAdiabaticFrameHamiltonian)(s::Real)
    res = zeros(eltype(h.h.m[1]), size(h.h.m[1]))
    update!(res, s, h.dθ)
    update!(res, s, h.h)
    res
end

function (h::UnitAdiabaticFrameHamiltonian)(t::Real)
    s = t / h._tf_cache[1]
    res = zeros(eltype(h.h.m[1]), size(h.h.m[1]))
    update!(res, s, h.dθ)
    update!(res, s, h.h)
    res
end

struct UnitlessAdiabaticFramePausingHamiltonian
    dθ::LinearOperator
    h::LinearOperator
    sp::Float64
    sd::Float64
    tf_ext::Float64
    n_qubit::Int64
end

struct UnitAdiabaticFramePausingHamiltonian
    dθ::LinearOperator
    h::LinearOperator
    sp::Float64
    sd::Float64
    tf::Float64
    tf_ext::Float64
    n_qubit::Int64
end


function construct_pausing_hamiltonian(sp, sd, H::UnitlessAdiabaticFrameHamiltonian)
    tf_ext = (1+sd)*H._tf_cache[1]
    h_op = multiply(1+sd, H.h)
    UnitlessAdiabaticFramePausingHamiltonian(H.dθ, h_op, sp, sd, tf_ext, H.n_qubit)
end

function construct_pausing_hamiltonian(sp, sd, H::UnitAdiabaticFrameHamiltonian)
    tf_ext = (1+sd)*H._tf_cache[1]
    dθ_op = multiply(1/(1+sd), H.dθ)
    UnitAdiabaticFramePausingHamiltonian(dθ_op, H.h, sp, sd, H._tf_cache[1], tf_ext, H.n_qubit)
end


function (h::UnitlessAdiabaticFramePausingHamiltonian)(s::Real)
    res = zeros(eltype(h.h.m[1]), size(h.h.m[1]))
    if s <= h.sp
        update!(res, s, h.dθ)
        update!(res, s, h.h)
    elseif s <= h.sp + h.sd
        update!(res, h.sp, h.h)
    else
        sn = s-h.sd
        update!(res, sn, h.dθ)
        update!(res, sn, h.h)
    end
    res
end

function (h::UnitAdiabaticFramePausingHamiltonian)(t::Real)
    s = t / h.tf
    res = zeros(eltype(h.h.m[1]), size(h.h.m[1]))
    if s <= h.sp
        update!(res, s, h.dθ)
        update!(res, s, h.h)
    elseif s <= h.sp + h.sd
        update!(res, h.sp, h.h)
    else
        sn = s-h.sd
        update!(res, sn, h.dθ)
        update!(res, sn, h.h)
    end
    res
end
