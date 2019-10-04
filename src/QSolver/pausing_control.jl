"""
$(TYPEDEF)
Pausing controller for annealing with adiabatic frame Hamiltonians.
# Fields
$(FIELDS)
"""
struct PausingControl <: AbstractAnnealingControl
    """Denotes times that the timestepping algorithm must step to. They should be where the pausing is turned on/off"""
    tstops
    """An function to convert the unitless time to annealing parameter s"""
    annealing_parameter
    PausingControl(t, a) = new(sort(t), a)
end


"""
    function (p::PausingControl)(tf::Real, t::Real)

Generate annealing parameter ``s``, scaling of adiabatic part and scaling of geometric part based on the type of `tf`. If `tf` is `Real`, the adiabatic part is scaled by `tf`.
"""
function (p::PausingControl)(tf::Real, t::Real)
    s = p.annealing_parameter(t)
    s, tf, 1.0
end


"""
    function (p::PausingControl)(tf::UnitTime, t::Real)

If `tf` is UnitTime object, the geometric part is scaled by `1/tf`.
"""
function (p::PausingControl)(tf::UnitTime, t::Real)
    s = p.annealing_parameter(t / tf)
    s, 1.0, 1 / tf
end


"""
    function adjust_u0(u0, p)

Convert the state vector/density matrix to the corresponding DEDataArray depending on the type of control `p`.
"""
function adjust_u0(u0, p::PausingControl)
    if ndims(u0) == 1
        DEPausingVec(u0, 1)
    elseif ndims(u0) == 2
        DEPausingMat(u0, 1)
    else
        throw(ArgumentError("u0 can either be a vector or matrix."))
    end
end


"""
$(TYPEDEF)
DEDataVector type used for pausing control.

# Fields
$(FIELDS)
"""
mutable struct DEPausingVec{T} <: DEDataVector{T}
    """The state vector"""
    x::Array{T,1}
    """Flag for whether pause is turned on"""
    pause::Int
end


"""
$(TYPEDEF)
DEDataMatrix type used for pausing control.

# Fields
$(FIELDS)
"""
mutable struct DEPausingMat{T} <: DEDataMatrix{T}
    """The density matrix"""
    x::Array{T,2}
    """Flag for whether pause is turned on"""
    pause::Int
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
    res = Array{eltype(tstops),1}()
    for t in tstops
        idx = findlast((x) -> x < t, starts)
        if idx == nothing
            push!(res, t)
        else
            push!(res, t + cum_inter[idx])
        end
    end
    append!(res, p.tstops)
    unique(sort(res))
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
    du.x .+= -1.0im * geometric_scale * u.pause * G * u.x
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
    du.x .+= -1.0im * geometric_scale * u.pause * (G * u.x - u.x * G)
end


function (h::AdiabaticFrameHamiltonian)(
    u::Union{DEPausingVec,DEPausingMat},
    p::AbstractAnnealingParams,
    t::Real
)
    s, adiabatic_scale, geometric_scale = p.control(p.tf, t)
    ω = 2.0 * π * adiabatic_scale * h.diagonal(s)
    G = geometric_scale * u.pause * h.geometric(s)
    ω + G
end


function (h::AdiabaticFrameHamiltonian)(
    u::Union{DEPausingVec,DEPausingMat},
    a_scale::Real,
    g_scale::Real,
    s::Real
)
    ω = 2.0 * π * a_scale * h.diagonal(s)
    G = g_scale * u.pause * h.geometric(s)
    ω + G
end


function single_pausing(sp, sd)
    function res(s)
        if s <= sp
            s
        elseif sp < s <= sp + sd
            sp
        elseif sp + sd < s
            s - sd
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


function prepare_callback(kw_dict, control::PausingControl)
    cb = DiscreteCallback(pause_condition, pause_affect!)
    res = Dict{Symbol,Any}(kw_dict)
    res[:callback] = cb
end


struct AdjustedTimeDependentCoupling
    coupling::TimeDependentCoupling
    annealing_parameter
end


function (c::AdjustedTimeDependentCoupling)(t)
    s = c.annealing_parameter(t)
    c.coupling(s)
end


function attach_annealing_param(p::PausingControl, c::TimeDependentCouplings)
    res = [attach_annealing_param(p, x) for x in c.coupling]
    TimeDependentCouplings(res...)
end


function attach_annealing_param(p::PausingControl, c::TimeDependentCoupling)
    AdjustedTimeDependentCoupling(c, p.annealing_parameter)
end


function (D::AFRWADiffEqOperator{S})(du, u, p, t) where S <: PausingControl
    s, a_scale, g_scale = p.control(p.tf, t)
    w, v = ω_matrix_RWA(D.H, u, p.tf, s, D.lvl)
    ρ = v' * u.x * v
    H = Diagonal(w)
    cache = -1.0im * a_scale * (H * ρ - ρ * H)
    ω_ba = repeat(w, 1, length(w))
    ω_ba = transpose(ω_ba) - ω_ba
    D.Davies(cache, ρ, ω_ba, v, p.tf, s)
    mul!(du, v, cache * v')
end


function ω_matrix_RWA(H::AdiabaticFrameHamiltonian, u::DEPausingMat, tf, s, lvl)
    ω = 2π * H.diagonal(s)
    off = 2π * u.pause * H.geometric(s) / tf
    ω + off
    eigen!(Hermitian(ω + off), 1:lvl)
end
