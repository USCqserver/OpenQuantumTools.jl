struct PausingControl <: AbstractAnnealingControl
    tstops
    annealing_parameter
end


function (p::PausingControl)(tf::Real, t::Real)
    s = p.annealing_parameter(t)
    s, tf, 1.0
end


function (p::PausingControl)(tf::UnitTime, t::Real)
    s = p.annealing_parameter(t / tf)
    s, 1.0, 1 / tf
end


function QTBase.adjust_sspan(p::PausingControl, sspan)
    starts = p.tstops[1:2:end]
    ends = p.tstops[2:2:end]
    sf = 1 + sum(ends - starts)
    (sspan[1] * sf, sspan[2] * sf)
end


function QTBase.adjust_tstops(p::PausingControl, tstops)
    pstop = sort(p.tstops)
    starts = pstop[1:2:end]
    ends = pstop[2:2:end]
    cum_inter = cumsum(ends .- starts)
    res = Array{eltype(tstops), 1}()
    for t in tstops
        idx = findlast((x)->x<t, starts)
        if idx == nothing
            push!(res, t)
        else
            push!(res, t+cum_inter[idx])
        end
    end
    append!(res, p.tstops)
    unique(sort(res))
end


mutable struct DEPausingVec{T} <: DEDataVector{T}
    x::Array{T,1}
    pause::Int
end


mutable struct DEPausingMat{T} <: DEDataMatrix{T}
    x::Array{T,2}
    pause::Int
end


function (h::AdiabaticFrameHamiltonian)(
    du::DEPausingVec,
    u::DEPausingVec,
    p::AbstractAnnealingParams,
    t::Real
)
    s, adiabatic_scale, geometric_scale = p.control(p.tf, t)
    ω = h.diagonal(s)
    du.x .= -2.0im * π * adiabatic_scale * ω * u.x
    G = h.geometric(s)
    du.x .+= -2.0im * π * geometric_scale * u.pause * G * u.x
end


function (h::AdiabaticFrameHamiltonian)(
    du::DEPausingMat,
    u::DEPausingMat,
    p::AbstractAnnealingParams,
    t::Real
)
    s, adiabatic_scale, geometric_scale = p.control(p.tf, t)
    ω = h.diagonal(s)
    du.x .= -2.0im * π * adiabatic_scale * (ω * u.x - u.x * ω)
    G = h.geometric(s)
    du.x .+= -2.0im * π * geometric_scale * u.pause * (G * u.x - u.x * G)
end


function (h::AdiabaticFrameHamiltonian)(
    u::Union{DEPausingVec, DEPausingMat},
    p::AbstractAnnealingParams,
    t::Real
)
    s, adiabatic_scale, geometric_scale = p.control(p.tf, t)
    ω = 2.0 * π * adiabatic_scale * h.diagonal(s)
    G = 2.0 * π * geometric_scale * u.pause * h.geometric(s)
    ω + G
end


function single_pausing(sp, sd)
    function res(s)
        if s <= sp
            s
        elseif sp < s <= sp + sd
            sp
        elseif sp + sd < s
            s - sp
        end
    end
    PausingControl([sp, sp + sd], res)
end


function pause_condition(u, t, integrator)
    if typeof(integrator.p.tf) <: Real
        return t in integrator.p.control.tstops
    elseif typeof(integrator.p.tf) <: UnitTime
        return t in integrator.p.tf * integrator.p.control.tstops
    end
end


function pause_affect!(integrator)
    for c in full_cache(integrator)
        c.pause = 1 - c.pause
    end
    u_modified!(integrator, false)
end
