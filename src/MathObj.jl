export σx, σz, σy, σi, σ, ⊗, PauliVec, comm, comm!, Hamiltonian
export matrix_decompose, check_positivity
export eigen_value_eval, eigen_state_eval, inst_population, gibbs_state

const σx = [0.0+0.0im 1; 1 0]
const σy = [0.0+0.0im -1.0im; 1.0im 0]
const σz = [1.0+0.0im 0; 0 -1]
const σi = [1.0+0.0im 0; 0 1]
const σ = [σx, σy, σz, σi]
xvec = [[1.0+0.0im, 1.0]/sqrt(2), [1.0+0.0im, -1.0]/sqrt(2)]
yvec = [[1.0im, -1.0]/sqrt(2), [1.0im, 1.0]/sqrt(2)]
zvec = [[1.0+0.0im, 0], [0, 1.0+0.0im]]
const PauliVec = [xvec, yvec, zvec]

⊗ = kron

"""
    comm(A, B)

Calculate the commutator of matrices `A` and `B`
"""
function comm(A, B)
    A*B - B*A
end

"""
    comm!(Y, A, B)

Calculate the commutator of matrices `A` and `B`. Write the result in `Y`
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

" Evaluate the Hamiltonian at time `t` "
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
    for i in range(2,stop=length(t))
        for j in range(1,length=res_dim)
            if real(res[i][j]'*res[i-1][j]) < 0
                res[i][j] = -res[i][j]
            end
        end
    end
    res
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
    matrix_decompose(mat::Array{T,2}, basis::Array{Array{T,2},1})

Decompse matrix `mat` onto matrix basis `basis`

# Examples
```julia-repl
julia> matrix_decompose(1.0*σx+2.0*σy+3.0*σz, [σx,σy,σz])
3-element Array{Complex{Float64},1}:
 1.0 + 0.0im
 2.0 + 0.0im
 3.0 + 0.0im
```
"""
function matrix_decompose(mat::Array{T,2}, basis::Array{Array{T,2},1}) where T<:Number
    dim = size(basis[1])[1]
    [tr(mat*b)/dim for b in basis]
end

"""
    check_positivity(m)

Check if matrix `m` is positive. Internally it compares the minimum eigen value of `m` to 0.
"""
function check_positivity(m::Array{T,2}) where T<:Number
    if !ishermitian(m)
        @warn "Input fails the numerical test for Hermitian matrix. Use the upper triangle to construct a new Hermitian matrix."
        d = Hermitian(m)
    else
        d = m
    end
    eigmin(d) > 0
end

"""
    gibbs_state(h, β)

Calculate the Gibbs state of the matrix `h` at temperature `β`.

# Examples
```julia-repl
julia> gibbs_state(σz, 0.1)
2×2 Array{Complex{Float64},2}:
 0.119203+0.0im       0.0+0.0im
      0.0+0.0im  0.880797+0.0im
```
"""
function gibbs_state(h, β)
    res = zeros(eltype(h), size(h))
    Z = 0.0
    eig_sys = eigen(Hermitian(h))
    for (i, E) in enumerate(eig_sys.values)
        t = exp(-β*E)
        BLAS.syr!('U', Complex(t), eig_sys.vectors[:, i], res)
        Z += t
    end
    res/Z
end
