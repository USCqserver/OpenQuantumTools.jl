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
