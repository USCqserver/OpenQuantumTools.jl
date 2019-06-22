"""
    LowLevelSystem

Object for a projected low level system. The projection is only valid for real Hamiltonians.

**Fields**
- `s` -- unitless time grid.
- `ev` -- energy values for different levels.
- `dθ` -- geometric terms.
- `op` -- projected system bath interaction operators
- `ref` -- energy eigenstates at the final time
- `lvl` -- number of levels being kept.
"""
struct LowLevelSystem
    s::AbstractArray{Float64, 1}
    ev::Array{Array{Float64, 1}, 1}
    dθ::Array{Array{Float64, 1}, 1}
    op::Array{Array{Array{Float64, 2},1},1}
    ref::Array{Float64, 2}
    lvl::Int
end

"""
    RotatedTwoLevelSystem

Object for a rotated two level system.

**Fields**
- `s` -- unitless time grid.
- `ω` -- ω₁₂.
- `T` -- off diagonal element of the Hamiltonian.
- `G` -- geometric term.
- `a`, `b`, `c`, `d` -- projected system bath parameters
- `θ` -- rotation angle
"""
struct RotatedTwoLevelSystem
    s::AbstractArray{Float64, 1}
    ω::Array{Float64, 1}
    T::Array{Float64, 1}
    G::Array{Float64, 1}
    a::Array{Float64, 1}
    b::Array{Float64, 1}
    c::Array{Float64, 1}
    d::Array{Float64, 1}
    θ::Array{Float64, 1}
end

"""
    get_dθ(sys::LowLevelSystem, i=1, j=2)

Get the geometric terms between i, j energy levels from `LowLevelSystem`.
"""
function get_dθ(sys::LowLevelSystem, i=1, j=2)
    if j > i
        idx = (2 * sys.lvl - i ) * (i - 1) ÷ 2 + (j - i)
        return [x[idx] for x in sys.dθ]
    elseif j < i
        idx = (2 * sys.lvl - j ) * (j - 1) ÷ 2 + (i - j)
        return [-x[idx] for x in sys.dθ]
    else
        error("No diagonal element for dθ.")
    end
end

"""
    proj_low_lvl(hfun, dhfun, interaction, s_axis::AbstractArray{T, 1}; ref=nothing, tol=1e-4, lvl=2) where T<:Number

Project a multi-level quantum system to the lowest `lvl` levels. `hfun` and `dhfun` are the Hamiltonian and its derivative. `interaction` is a list of system bath interaction operators. `s_axis` is the grid for the computation. `ref` is the reference for the initial eigenstates. `tol` specify the error tolerance for eigen-decomposition if sparse matrix is used.
"""
function proj_low_lvl(hfun, dhfun, interaction, s_axis::AbstractArray{T, 1}; ref=nothing, tol=1e-4, lvl=2) where T<:Number
    H = hfun(s_axis[1])
    dH = dhfun(s_axis[1])
    ev = Array{Float64, 1}()
    dθ = Array{Array{Float64, 1}, 1}()
    op = Array{Array{Array{Float64, 2},1},1}()
    if ref==nothing
        low_obj = LowLevelSystem(s_axis, ev, dθ, op, zeros(size(H,1), lvl), lvl)
        # v0 for eigs can not be a full zero vector, need to deal with separately
        if issparse(H)
            w, v = eigs(H, nev=lvl, which=:SR, tol=tol)
        else
            w, v = eigen!(H)
            w = w[1:lvl]
            v = v[:, 1:lvl]
        end
        params_push!(low_obj, w, v, dH, interaction)
    else
        low_obj = LowLevelSystem(s_axis, ev, dθ, op, ref, size(ref, 2))
        _proj_step!(low_obj, H, dH, interaction, lvl, tol)
    end

    for s in s_axis[2:end]
        H = hfun(s)
        dH = dhfun(s)
        _proj_step!(low_obj, H, dH, interaction, lvl, tol)
    end
    low_obj
end

function _proj_step!(params::LowLevelSystem, H::SparseMatrixCSC{T, V}, dH, interaction, lvl, tol) where T<:Number where V<:Int
    w, v = eigs(H, nev=lvl, which=:SR, tol=tol, v0=params.ref[:,1])
    params_push!(params, w, v, dH, interaction)
end

function _proj_step!(params::LowLevelSystem, H::Array{T, 2}, dH, interaction, lvl, tol) where T<:Number
    w, v = eigen!(H)
    params_push!(params, w[1:lvl], v[:, 1:lvl], dH, interaction)
end

function params_push!(params::LowLevelSystem, w, v, dH, interaction)
    push!(params.ev, w)

    # update reference vectors
    for i in 1:params.lvl
        if v[:, i]' * params.ref[:, i] < 0
            params.ref[:, i] = -v[:, i]
        else
            params.ref[:, i] = v[:, i]
        end
    end

    # update dθ
    res = []
    for i in 1:params.lvl
        for j in 1+i:params.lvl
            vi = @view params.ref[:, i]
            vj = @view params.ref[:, j]
            t = vi' * dH * vj / (w[i] - w[j])
            push!(res, t)
        end
    end
    push!(params.dθ, res)

    # update projected interaction operators
    op = [params.ref'*x*params.ref for x in interaction]
    push!(params.op, op)
end
