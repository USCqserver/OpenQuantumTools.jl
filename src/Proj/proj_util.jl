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
    get_dθ(sys::LowLevelSystem[, i=1, j=2])

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


function _push_dθ!(params::LowLevelSystem, dH, w)
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
end

function _update_ref!(params::LowLevelSystem, v)
    for i in 1:params.lvl
        if v[:, i]' * params.ref[:, i] < 0
            params.ref[:, i] = -v[:, i]
        else
            params.ref[:, i] = v[:, i]
        end
    end
end

function _init_lowlevel_params(s_axis, w, v, dH, interaction, lvl)
    l = length(s_axis)
    ev = Array{Array{Float64, 1}, 1}()
    push!(ev, w)
    dθ = Array{Array{Float64, 1}, 1}()
    dθt = []
    for i in 1:lvl
        for j in 1+i:lvl
            vi = @view v[:, i]
            vj = @view v[:, j]
            t = vi' * dH * vj / (w[i] - w[j])
            push!(dθt, t)
        end
    end
    push!(dθ, dθt)
    op = Array{Array{Array{Float64, 2},1},1}()
    opt = [v'*x*v for x in interaction]
    push!(op, opt)
    ref = v
    LowLevelSystem(s_axis, ev, dθ, op, ref, lvl)
end

function empty_lowlevel_params(s_axis, ref)
    ev = Array{Float64, 1}()
    dθ = Array{Array{Float64, 1}, 1}()
    op = Array{Array{Array{Float64, 2},1},1}()
    LowLevelSystem(s_axis, ev, dθ, op, ref, size(ref, 2))
end

function params_push!(params::LowLevelSystem, w, v, dH, interaction)
    push!(params.ev, w)
    _update_ref!(params, v)
    _push_dθ!(params, dH, w)
    op = [params.ref'*x*params.ref for x in interaction]
    push!(params.op, op)
end

function _init_proj_step(H::SparseMatrixCSC{T, V}, dH, interaction, s_axis, lvl, tol) where T<:Number where V<:Int
    w, v = eigs(H, nev=lvl, which=:SR, tol=tol)
    _init_lowlevel_params(s_axis, w, v, dH, interaction, lvl)
end

function _init_proj_step(H::Array{T, 2}, dH, interaction, s_axis, lvl, tol) where T<:Number
    w, v = eigen!(H)
    _init_lowlevel_params(s_axis, w[1:lvl], v[:, 1:lvl], dH, interaction, lvl)
end

function _proj_step!(params::LowLevelSystem, H::SparseMatrixCSC{T, V}, dH, interaction, lvl, tol) where T<:Number where V<:Int
    w, v = eigs(H, nev=lvl, which=:SR, tol=tol, v0=params.ref[:,1])
    params_push!(params, w, v, dH, interaction)
end

function _proj_step!(params::LowLevelSystem, H::Array{T, 2}, dH, interaction, lvl, tol) where T<:Number
    w, v = eigen!(H)
    params_push!(params, w[1:lvl], v[:, 1:lvl], dH, interaction)
end

function proj_low_lvl(hfun, dhfun, interaction, s_axis::AbstractArray{T, 1}; ref=nothing, tol=1e-4, lvl=2) where T<:Number
    H = hfun(s_axis[1])
    dH = dhfun(s_axis[1])
    if ref==nothing
        low_obj = _init_proj_step(H, dH, interaction, s_axis, lvl, tol)
    else
        low_obj = empty_lowlevel_params(s_axis, ref)
        _proj_step!(low_obj, H, dH, interaction, lvl, tol)
    end
    for s in s_axis[2:end]
        H = hfun(s)
        dH = dhfun(s)
        _proj_step!(low_obj, H, dH, interaction, lvl, tol)
    end
    low_obj
end

function optimal_interaction_angle(low::LowLevelSystem)
    if low.lvl!=2
        error("Optimal rotation only works for the lowest 2 levels.")
    end

    function obj(θ, op)
        g = [1.0, 0]
        e = [0, 1.0]
        pointer_g = -cos(θ/2 - π/4) * g + cos(θ/2 + π/4) * e
        pointer_e = -cos(θ/2 + π/4) * g - cos(θ/2 - π/4) * e
        res = 0.0
        for o in op
            res += abs2(pointer_e' * o * pointer_e - pointer_g' * o * pointer_g)
        end
        -res
    end

    opt_θ = Array{Float64, 1}()
    for op in low.op
        f = (x)->obj(x, op)
        opt = optimize(f, -π/2, π/2)
        push!(opt_θ, opt.minimizer)
    end
    opt_θ
end

function _rotate_by_interaction(sys::LowLevelSystem)
    θ = optimal_interaction_angle(sys)
    g = [1.0, 0]
    e = [0, 1.0]
    ω = []
    T = []
    a = []
    b = []
    c = []
    d = []
    for i in eachindex(sys.s)
        H = Diagonal([sys.ev[i][1], sys.ev[i][2]])
        U = [-cos(θ[i]/2 - π/4) -cos(θ[i]/2 + π/4); cos(θ[i]/2 + π/4) -cos(θ[i]/2 - π/4)]
        Hr = U' * H * U
        push!(ω, real(Hr[1,1]-Hr[2,2]))
        push!(T, Hr[1,2])
        at = 0.0
        bt = 0.0
        ct = 0.0
        dt = 0.0
        for op in sys.op[i]
            or = U' * op * U
            at += (or[1,1] - or[2,2])^2
            bt += abs2(or[1,2])
            ct += or[1,2] * (or[1,1] - or[2,2])
            dt +=  or[1,2] * (or[1,1] + or[2,2])
        end
        push!(a, at)
        push!(b, bt)
        push!(c, ct)
        push!(d, dt)
    end
    θ_itp = construct_interpolations(sys.s, θ, extrapolation="line")
    θᴱ = gradient(θ_itp, sys.s) / 2
    dθ = get_dθ(sys, 1, 2)
    G = dθ - θᴱ
    RotatedTwoLevelSystem(sys.s, ω, T, G, a, b, c, d, θ)
end

function _rotate_by_LZ(sys::LowLevelSystem)
    dθ_itp = construct_interpolations(sys.s, get_dθ(sys), extrapolation="line")
    θᴸ = [quadgk(dθ_itp, 0, s)[1] for s in sys.s]
    g = [1.0, 0]
    e = [0, 1.0]
    ω = []
    T = []
    a = []
    b = []
    c = []
    d = []
    for i in eachindex(sys.s)
        H = Diagonal([sys.ev[i][1], sys.ev[i][2]])
        U = [cos(θᴸ[i]) -sin(θᴸ[i]); sin(θᴸ[i]) cos(θᴸ[i])]
        Hr = U' * H * U
        push!(ω, real(Hr[1,1]-Hr[2,2]))
        push!(T, Hr[1,2])
        at = 0.0
        bt = 0.0
        ct = 0.0
        dt = 0.0
        for op in sys.op[i]
            or = U' * op * U
            at += (or[1,1] - or[2,2])^2
            bt += abs2(or[1,2])
            ct += or[1,2] * (or[1,1] - or[2,2])
            dt +=  or[1,2] * (or[1,1] + or[2,2])
        end
        push!(a, at)
        push!(b, bt)
        push!(c, ct)
        push!(d, dt)
    end
    G = zeros(length(sys.s))
    RotatedTwoLevelSystem(sys.s, ω, T, G, a, b, c, d, θᴸ)
end

function _zero_rotate(sys::LowLevelSystem)
    g = [1.0, 0]
    e = [0, 1.0]
    ω = []
    G = []
    a = []
    b = []
    c = []
    d = []
    for i in eachindex(sys.s)
        push!(ω, sys.ev[i][1]-sys.ev[i][2])
        push!(G, sys.dθ[i][1])
        at = 0.0
        bt = 0.0
        ct = 0.0
        dt = 0.0
        for op in sys.op[i]
            at += (op[1,1] - op[2,2])^2
            bt += abs2(op[1,2])
            ct += op[1,2] * (op[1,1] - op[2,2])
            dt +=  op[1,2] * (op[1,1] + op[2,2])
        end
        push!(a, at)
        push!(b, bt)
        push!(c, ct)
        push!(d, dt)
    end
    T = zeros(length(sys.s))
    RotatedTwoLevelSystem(sys.s, ω, T, G, a, b, c, d, zeros((0,)))
end

"""
    rotate_lowlevel_system(sys::LowLevelSystem; method=nothing)

Rotate the projected 2 level system in adiabatic frame to a different rotating frame according to `method`.

**method**
- `"none"` (default) -- no ratation at all.
- `"env"` -- rotate the system to maximize distinguishability of states with respect to system bath coupling.
- `lz` -- rotate the system to minimize the Landau-Zener transistion.
"""
function rotate_lowlevel_system(sys::LowLevelSystem; method="none")
    if method == "none"
        return _zero_rotate(sys)
    elseif method == "env"
        return _rotate_by_interaction(sys)
    elseif method == "lz"
        return _rotate_by_LZ(sys)
    else
        @warn "No specific method: " method
    end
end
