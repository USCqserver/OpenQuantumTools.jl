"""
    rotate_lowlevel_system(sys::LowLevelSystem; method=nothing)

Rotate the projected 2 level system in adiabatic frame to a different rotating frame according to `method`.

**method**
- `"none"` (default) -- no ratation at all.
- `"env"` -- rotate the system to maximize distinguishability of states with respect to system bath coupling.
- `lz` -- rotate the system to minimize the Landau-Zener transistion.
"""
function rotate_lowlevel_system(sys::LowLevelSystem; method="none", rotation_point=nothing)
    if method == "none"
        return _zero_rotate(sys)
    elseif method == "env"
        return _rotate_by_interaction(sys)
    elseif method == "lz"
        return _rotate_by_LZ(sys)
    elseif method == "composite"
        return _composite_rotate(sys, rotation_point)
    else
        @warn "No specific method: " method
    end
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

function optimal_interaction_angle(low::LowLevelSystem)
    if low.lvl!=2
        error("Optimal rotation only works for the lowest 2 levels.")
    end

    # objection function for optimization
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

function _composite_rotate(sys::LowLevelSystem, rotation_point)
    dθ_itp = construct_interpolations(sys.s[rotation_point:end], get_dθ(sys)[rotation_point:end], extrapolation="line")
    θᴸ_2 = [quadgk(dθ_itp, sys.s[rotation_point], s)[1] for s in sys.s[rotation_point:end]]
    θᴸ_1 = zeros(rotation_point-1)
    θᴸ = vcat(θᴸ_1, θᴸ_2)
    g = [1.0, 0]
    e = [0, 1.0]
    ω = []
    T = []
    G = []
    a = []
    b = []
    c = []
    d = []
    # === first half (no rotation) ===
    for i in 1:rotation_point-1
        push!(ω, sys.ev[i][1]-sys.ev[i][2])
        push!(G, sys.dθ[i][1])
        push!(T, 0.0)
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
    # === second half ===
    for i in rotation_point:length(sys.s)
        H = Diagonal([sys.ev[i][1], sys.ev[i][2]])
        U = [cos(θᴸ[i]) -sin(θᴸ[i]); sin(θᴸ[i]) cos(θᴸ[i])]
        Hr = U' * H * U
        push!(ω, real(Hr[1,1]-Hr[2,2]))
        push!(T, Hr[1,2])
        push!(G, 0.0)
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
    RotatedTwoLevelSystem(sys.s, ω, T, G, a, b, c, d, θᴸ)
end
