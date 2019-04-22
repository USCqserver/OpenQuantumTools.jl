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
    comm(A, B)

Calculate the commutator of matrices `A` and `B`.
"""
function comm(A, B)
    A*B - B*A
end

"""
    comm!(Y, A, B)

Calculate the commutator of matrices `A` and `B` and add the result to `Y`.
"""
function comm!(Y, A, B)
    mul!(Y,A,B)
    axpy!(-1.0,B*A,Y)
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
    eigen_state_continuation!(eigen_states, reference=nothing)

Give an `eigen_states` matrix with dimension (dim, dim, time), adjust the sign of each eigenstate such that the inner product between the neighboring eigenstates are positive.
"""
function eigen_state_continuation!(eigen_states, reference=nothing)
    dim = size(eigen_states)
    if reference != nothing
        for i in range(1, length=dim[2])
            if real(eigen_states[:, i, 1]' * reference[:, i]) < 0
                eigen_states[:, i, 1] = -eigen_states[:, i, 1]
            end
        end
    end
    for i in range(2, stop=size[3])
        for j in range(1, length=dim[2])
            if real(eigen_states[:, j, i]' * eigen_states[:, j, i-1]) < 0
                eigen_states[:, j, i] = -eigen_states[:, j, i]
            end
        end
    end
end


"""
    eigen_eval(hfun, t; levels=2, tol=1e-4)

Calculate the eigen values and eigen states of Hamiltonian `hfun` at each points of `t`. The output keeps the lowest `levels` eigenstates and their corresponding eigenvalues. `tol` specifies the error tolerance for sparse matrices decomposition.
"""
function eigen_eval(hfun, t::AbstractArray{Float64,1}; levels::Int=2, tol=1e-4)
    t_dim = length(t)
    H = hfun(t[1])
    res_val = Array{eltype(H), 2}(undef, (levels, t_dim))
    res_vec = Array{eltype(H), 3}(undef, (size(H)[1], levels, t_dim))
    if issparse(H)
        eigfun = (x)-> eigs(x, nev=levels, which=:SR, tol=tol)
    else
        eigfun = (x)-> eigen(Hermitian(x))
    end
    val, vec = eigfun(H)
    res_val[:, 1] = val[1:levels]
    res_vec[:, :, 1] = vec[:, 1:levels]
    for (i, t_val) in enumerate(t[2:end])
        val, vec = eigfun(hfun(t_val))
        res_val[:, i+1] = val[1:levels]
        res_vec[:, :, i+1] = vec[:, 1:levels]
    end
    res_val, res_vec
end

function proj_2lvl(hfun, dh, t_list)

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
```
"""
function gibbs_state(h, β)
    res = zeros(ComplexF64, size(h))
    Z = 0.0
    eig_sys = eigen(Hermitian(h))
    for (i, E) in enumerate(eig_sys.values)
        t = exp(-β*E)
        her!('U', t, eig_sys.vectors[:, i], res)
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
```
"""
function low_level_hamiltonian(h, levels)
    if levels > size(h)[1]
        @warn "Subspace dimension bigger than total dimension."
        return h
    else
        eig_sys = eigen(Hermitian(h))
        res = zeros(eltype(h), size(h))
        for i in range(1,stop=levels)
            her!('U', eig_sys.values[i], eig_sys.vectors[:, i], res)
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

function proj_2lvl(hfun, dhfun, interaction, s_axis::AbstractArray{T, 1}; reference=nothing, tol=1e-4) where T<:Number

    ev = Array{Float64, 2}(undef, (2, length(s_axis)))
    dθ = Array{Float64, 1}(undef, length(s_axis))
    op = Array{Float64, 2}(undef, (3, length(s_axis)))

    H = hfun(s_axis[1])
    dH = dhfun(s_axis[1])

    if issparse(H)
        # sparse case
        if reference==nothing
            w, pv = eigs(H, nev=2, which=:SR, tol=tol)
        else
            w, pv = eigs(H, nev=2, which=:SR, tol=tol, v0=reference[:,1])
            for j in 1:2
                if pv[:, j]' * reference[:, j] < 0
                    pv[:, j] = -pv[:, j]
                end
            end
        end
        ev[:, 1] = w
        p1 = @view pv[:,1]
        p2 = @view pv[:,2]
        dθ[1] = p2' * dH * p1 / (w[2] - w[1])
        op[1, 1] = p1' * interaction * p1
        op[2, 1] = p2' * interaction * p2
        op[3, 1] = p2' * interaction * p1

        for (i, s) in enumerate(s_axis[2:end])
            H = hfun(s)
            dH = dhfun(s)
            w, cv = eigs(H, nev=2, which=:SR, tol=tol, v0=pv[:, 1])
            for j in 1:2
                if cv[:, j]' * pv[:, j] < 0
                    pv[:, j] = -cv[:, j]
                else
                    pv[:, j] = cv[:, j]
                end
            end
            ev[:, i+1] = w
            p1 = @view pv[:,1]
            p2 = @view pv[:,2]
            dθ[i+1] = p2' * dH * p1 / (w[2] - w[1])
            op[1, i+1] = p1' * interaction * p1
            op[2, i+1] = p2' * interaction * p2
            op[3, i+1] = p2' * interaction * p1
        end

    else
        # dense case
        w, pv = eigen!(Hermitian(H))
        if reference!=nothing
            for j in 1:2
                if pv[:, j]' * reference[:, j] < 0
                    pv[:, j] = -pv[:, j]
                end
            end
        end
        ev[:, 1] = w
        p1 = @view pv[:,1]
        p2 = @view pv[:,2]
        dθ[1] = p2' * dH * p1 / (w[2] - w[1])
        op[1, 1] = p1' * interaction * p1
        op[2, 1] = p2' * interaction * p2
        op[3, 1] = p2' * interaction * p1

        for (i, s) in enumerate(s_axis[2:end])
            H = hfun(s)
            dH = dhfun(s)
            w, cv = eigen!(Hermitian(H))
            for j in 1:2
                if cv[:, j]' * pv[:, j] < 0
                    pv[:, j] = -cv[:, j]
                else
                    pv[:, j] = cv[:, j]
                end
            end
            ev[:, i+1] = w
            p1 = @view pv[:,1]
            p2 = @view pv[:,2]
            dθ[i+1] = p2' * dH * p1 / (w[2] - w[1])
            op[1, i+1] = p1' * interaction * p1
            op[2, i+1] = p2' * interaction * p2
            op[3, i+1] = p2' * interaction * p1
        end
    end
    ev, dθ, op, pv
end

function proj_2lvl_polaron(hfun, dhfun, interaction, s_axis::AbstractArray{T, 1}; reference=nothing, tol=1e-4) where T<:Number

    ev = Array{Float64, 2}(undef, (2, length(s_axis)))
    dθ = Array{Float64, 1}(undef, length(s_axis))
    op = Array{Float64, 2}(undef, (8, length(s_axis)))

    H = hfun(s_axis[1])
    dH = dhfun(s_axis[1])

    if issparse(H)
        # sparse case
        if reference==nothing
            w, pv = eigs(H, nev=2, which=:SR, tol=tol)
        else
            w, pv = eigs(H, nev=2, which=:SR, tol=tol, v0=reference[:,1])
            for j in 1:2
                if pv[:, j]' * reference[:, j] < 0
                    pv[:, j] = -pv[:, j]
                end
            end
        end
        ev[:, 1] = w
        p1 = @view pv[:,1]
        p2 = @view pv[:,2]
        dθ[1] = p2' * dH * p1 / (w[2] - w[1])
        op[:, 1] = _decompose_interaction_abcd(pv, interaction)

        for (i, s) in enumerate(s_axis[2:end])
            H = hfun(s)
            dH = dhfun(s)
            w, cv = eigs(H, nev=2, which=:SR, tol=tol, v0=pv[:, 1])
            for j in 1:2
                if cv[:, j]' * pv[:, j] < 0
                    pv[:, j] = -cv[:, j]
                else
                    pv[:, j] = cv[:, j]
                end
            end
            ev[:, i+1] = w
            p1 = @view pv[:,1]
            p2 = @view pv[:,2]
            dθ[i+1] = p2' * dH * p1 / (w[2] - w[1])
            op[:, i+1] = _decompose_interaction_abcd(pv, interaction)
        end

    else
        # dense case
        w, pv = eigen!(Hermitian(H))
        if reference!=nothing
            for j in 1:2
                if pv[:, j]' * reference[:, j] < 0
                    pv[:, j] = -pv[:, j]
                end
            end
        end
        ev[:, 1] = w
        p1 = @view pv[:,1]
        p2 = @view pv[:,2]
        dθ[1] = p2' * dH * p1 / (w[2] - w[1])
        op[:, 1] = _decompose_interaction_abcd(pv, interaction)

        for (i, s) in enumerate(s_axis[2:end])
            H = hfun(s)
            dH = dhfun(s)
            w, cv = eigen!(Hermitian(H))
            for j in 1:2
                if cv[:, j]' * pv[:, j] < 0
                    pv[:, j] = -cv[:, j]
                else
                    pv[:, j] = cv[:, j]
                end
            end
            ev[:, i+1] = w
            p1 = @view pv[:,1]
            p2 = @view pv[:,2]
            dθ[i+1] = p2' * dH * p1 / (w[2] - w[1])
            op[:, i+1] = _decompose_interaction_abcd(pv, interaction)
        end
    end
    ev, dθ, op, pv
end



function _decompose_interaction_abcd(v, inter_op)
    p1 = @view v[:, 1]
    p2 = @view v[:, 2]
    σ11 = Array{Float64, 1}()
    σ22 = Array{Float64, 1}()
    σ12 = Array{Float64, 1}()
    for op in inter_op
        push!(σ11, p1'*op*p1)
        push!(σ22, p2'*op*p2)
        push!(σ12, p1'*op*p2)
    end
    a12 = sum((σ11 - σ22).^2)
    c12 = sum(σ12 .* (σ11 - σ22))
    b11 = sum(abs2.(σ11))
    b22 = sum(abs2.(σ22))
    b12 = sum(abs2.(σ12))
    d11 = 2*sum(σ11 .* σ11)
    d22 = 2*sum(σ22 .* σ22)
    d12 = sum(σ12 .* (σ11 + σ22))
    [a12, c12, b11, b22, b12, d11, d22, d12]
end