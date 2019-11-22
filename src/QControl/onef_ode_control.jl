"""
$(TYPEDEF)

Defines a single fluctuator ensemble controller

# Fields

$(FIELDS)
"""
mutable struct FluctuatorControl{T} <: AbstractAnnealingControl
    """waitting time distribution for every fluctuators"""
    dist
    """cache for each fluctuator value"""
    b0
    """index of the fluctuator to be flipped next"""
    next_idx
    """time interval for next flip event"""
    next_τ
end


function FluctuatorControl(tf, num::Int, E::EnsembleFluctuator)
    dist = construct_distribution(tf, E)
    b0 = [x.b for x in E.f] .* rand([-1, 1], length(dist), num)
    next_τ, next_idx = findmin(rand(dist, num))
    FluctuatorControl{num}(dist, b0, next_idx, next_τ)
end


function (f::FluctuatorControl)()
    view(sum(f.b0, dims=1), :)
end


function next_state!(f::FluctuatorControl)
    next_τ, next_idx = findmin(rand(f.dist, size(f.b0, 2)))
    f.next_τ = next_τ
    f.next_idx = next_idx
    f.b0[next_idx] *= -1
    nothing
end


"""
$(TYPEDEF)

Defines stochastic system-bath coupling operator

# Fields

$(FIELDS)
"""
struct StochasticNoise <: AbstractOpenSys
    """system-bath coupling operator"""
    ops::AbstractCouplings
end


function (S::StochasticNoise)(A, u, tf::Real, t)
    A .+= tf * sum(u.n .* S.ops(t))
end


function (S::StochasticNoise)(A, u, tf::UnitTime, t)
    A .+= sum(u.n .* S.ops(t / tf))
end
