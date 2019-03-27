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

"""
    eigen_state_eval(hfun, t[, levels=[1,]])

Calculate the eigen states of Hamiltonian `hfun` at each points of `t`. The output levels are specified by `levels` argument.
"""
function eigen_state_eval(hfun, t::AbstractArray{Float64,1}; levels::Array{Int64,1}=[1,])
    res_dim = length(levels)
    res = Array{Array{Array{ComplexF64,1},1},1}(undef, length(t))
    for i in eachindex(t)
        eig_sys = eigen(Hermitian(hfun(t[i])))
        res[i] = [eig_sys.vectors[:,x] for x in levels]
    end
    res
end

"""
    eigen_state_continuation!(eigen_states)

Give a list of `eigen_states` at different time, adjust the sign of each eigen state such that the inner product between the neighboring eigen states are positive.
"""
function eigen_state_continuation!(eigen_states)
    for i in range(2,stop=length(eigen_states))
        for j in range(1,length=length(eigen_states[1]))
            if real(eigen_states[i][j]'*eigen_states[i-1][j]) < 0
                eigen_states[i][j] = -eigen_states[i][j]
            end
        end
    end
end

"""
    eigen_value_eval(hfun, t[, levels])

Calculate the eigen values of Hamiltonian `hfun` at each points of `t`. The output levels are specified by `levels` argument.
"""
function eigen_value_eval(hfun, t::AbstractArray{Float64,1}; levels::Array{Int64,1}=[1,])
    res_dim = length(levels)
    res = Array{Array{Float64,1},1}(undef, length(t))
    for i in eachindex(t)
        eig_sys = eigen(Hermitian(hfun(t[i])))
        res[i] = [eig_sys.values[x] for x in levels]
    end
    res
end

"""
    eigen_sys_eval(hfun, t[, levels])

Calculate the eigen values and eigen states of Hamiltonian `hfun` at each points of `t`. The output levels are specified by `levels` argument.
"""
function eigen_sys_eval(hfun, t::AbstractArray{Float64,1}; levels::Array{Int64,1}=[1,])
    res_dim = length(levels)
    res_value = Array{Array{Float64,1},1}(undef, length(t))
    res_vector = Array{Array{Array{ComplexF64,1},1},1}(undef, length(t))
    for i in eachindex(t)
        eig_sys = eigen(Hermitian(hfun(t[i])))
        res_value[i] = [eig_sys.values[x] for x in levels]
        res_vector[i] = [eig_sys.vectors[:,x] for x in levels]
    end
    res_value, res_vector
end

function sp_eigen_sys_eval(hfun, t::AbstractArray{Float64,1}; is_real=true, nev=2, which=:SR, tol=1e-4)
    if is_real == true
        res_vector = Array{Array{Float64, 2}, 1}(undef, length(t))
        res_value = Array{Array{Float64,1},1}(undef, length(t))
    else
        res_vector = Array{Array{ComplexF64, 2}, 1}(undef, length(t))
        res_value = Array{Array{ComplexF64,1},1}(undef, length(t))
    end
    for (idx, time) in enumerate(t)
        H =hfun(time)
        v, Φ = eigs(H, nev=nev, which=which, tol=tol)
        res_value[idx] = v
        res_vector[idx] = Φ
    end
    res_value, res_vector
end

"""
    function inst_population(t, states, hamiltonian; level=1)

For a time series quantum states given by `states`, whose time points are given by `t`, calculate the population of instantaneous eigenstates of `hamiltonian`. The levels of the instantaneous eigenstates are specified by `level`, which can be any slice index.
"""
function inst_population(t, states, hamiltonian; level=1:1)
    if typeof(level)<:Int
        level = level:level
    end
    pop = Array{Array{Float64, 1}, 1}(undef, length(t))
    for (i, v) in enumerate(t)
        hmat = hamiltonian(v)
        eig_sys = eigen(Hermitian(hmat))
        if ndims(states[i]) == 1
            inst_state = view(eig_sys.vectors,:,level)'
            pop[i] = abs2.(inst_state * states[i])
        elseif ndims(states[i]) == 2
            l = length(level)
            temp = Array{Float64, 1}(undef, l)
            for j in range(1, length=l)
                inst_state = view(eig_sys.vectors,:,j)
                temp[j] = real(inst_state'*states[i]*inst_state)
            end
            pop[i] = temp
        end
    end
    pop
end


"""
    gibbs_state(h, β)

Calculate the Gibbs state of the matrix `h` at temperature `β`.

# Examples
```julia-repl
julia> gibbs_state(σz, 0.1)
2×2 LinearAlgebra.Hermitian{Complex{Float64},Array{Complex{Float64},2}}:
 0.450166+0.0im       0.0+0.0im
      0.0-0.0im  0.549834+0.0im
"""
function gibbs_state(h, β)
    res = zeros(ComplexF64, size(h))
    Z = 0.0
    eig_sys = eigen(Hermitian(h))
    for (i, E) in enumerate(eig_sys.values)
        t = exp(-β*E)
        BLAS.her!('U', t, eig_sys.vectors[:, i], res)
        Z += t
    end
    Hermitian(res/Z)
end

"""
    low_level_hamiltonian(h, levels)

Calculate the Hamiltonian `h` projected to lower energy subspace containing `levels` energy levels.

# Examples
```julia-repl
julia> low_level_hamiltonian(σx⊗σx, 2)
4×4 LinearAlgebra.Hermitian{Complex{Float64},Array{Complex{Float64},2}}:
 -0.5+0.0im   0.0+0.0im   0.0+0.0im   0.5+0.0im
  0.0-0.0im  -0.5+0.0im   0.5+0.0im   0.0+0.0im
  0.0-0.0im   0.5-0.0im  -0.5+0.0im   0.0+0.0im
  0.5-0.0im   0.0-0.0im   0.0-0.0im  -0.5+0.0im
"""
function low_level_hamiltonian(h, levels)
    if levels > size(h)[1]
        @warn "Subspace dimension bigger than total dimension."
        return h
    else
        eig_sys = eigen(Hermitian(h))
        res = zeros(eltype(h), size(h))
        for i in range(1,stop=levels)
            BLAS.her!('U', eig_sys.values[i], eig_sys.vectors[:, i], res)
        end
        return Hermitian(res)
    end
end

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
