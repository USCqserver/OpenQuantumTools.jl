struct LowLevelParams
    s::AbstractArray{Float64, 1}
    ev::Array{Array{Float64, 1}, 1}
    dθ::Array{Array{Float64, 1}, 1}
    op::Array{Array{Array{Float64, 2},1},1}
    ref::Array{Float64, 2}
    lvl::Int
end

struct RotatedTwoLevelParams
    s::AbstractArray{Float64, 1}
    ω::Array{Float64, 1}
    T::Array{ComplexF64, 1}
    a::Array{Float64, 1}
    b::Array{Float64, 1}
    c::Array{Float64, 1}
    d::Array{Float64, 1}
end

function get_dθ(params::LowLevelParams, i, j)
    if j > i
        idx = (2 * params.lvl - i ) * (i - 1) ÷ 2 + (j - i)
        return [x[idx] for x in params.dθ]
    elseif j<i
        idx = (2 * params.lvl - j ) * (j - 1) ÷ 2 + (i - j)
        return [-x[idx] for x in params.dθ]
    else
        error("No diagonal element for dθ.")
    end
end

function _push_dθ!(params::LowLevelParams, dH, w)
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

function _update_ref!(params::LowLevelParams, v)
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
    LowLevelParams(s_axis, ev, dθ, op, ref, lvl)
end

function empty_lowlevel_params(s_axis, ref)
    ev = Array{Float64, 1}()
    dθ = Array{Array{Float64, 1}, 1}()
    op = Array{Array{Array{Float64, 2},1},1}()
    LowLevelParams(s_axis, ev, dθ, op, ref, size(ref, 2))
end

function params_push!(params::LowLevelParams, w, v, dH, interaction)
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

function _proj_step!(params::LowLevelParams, H::SparseMatrixCSC{T, V}, dH, interaction, lvl, tol) where T<:Number where V<:Int
    w, v = eigs(H, nev=lvl, which=:SR, tol=tol, v0=params.ref[:,1])
    params_push!(params, w, v, dH, interaction)
end

function _proj_step!(params::LowLevelParams, H::Array{T, 2}, dH, interaction, lvl, tol) where T<:Number
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

function optimal_interaction_rotation(low::LowLevelParams)
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

    opt_θ = []
    for op in low.op
        f = (x)->obj(x, op)
        opt = optimize(f, -π/2, π/2)
        push!(opt_θ, opt.minimizer)
    end
    opt_θ
end

function rotate_by_interaction(params::LowLevelParams, θ)
    g = [1.0, 0]
    e = [0, 1.0]
    ω = []
    T = []
    a = []
    b = []
    c = []
    d = []
    for i in eachindex(params.s)
        H = [params.ev[i][1] -1.0im*params.dθ[i][1]; 1.0im*params.dθ[i][1] params.ev[i][2]]
        U = [-cos(θ[i]/2 - π/4) -cos(θ[i]/2 + π/4); cos(θ[i]/2 + π/4) -cos(θ[i]/2 - π/4)]
        Hr = U * H * U'
        push!(ω, real(Hr[1,1]-Hr[2,2]))
        push!(T, Hr[1,2])
        at = 0.0
        bt = 0.0
        ct = 0.0
        dt = 0.0
        for op in params.op[i]
            or = U * op * U'
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
    RotatedTwoLevelParams(params.s, ω, T, a, b, c, d)
end
