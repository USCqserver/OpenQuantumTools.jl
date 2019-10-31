"""
$(TYPEDEF)

A symetric random telegraph noise with switch rate `γ/2` magnitude `b`

$(FIELDS)
"""
struct SymetricRTN
    "Magnitude"
    b
    "Two times the switching probability"
    γ
end


function correlation(τ, R::SymetricRTN)
    R.b^2*exp(-R.γ*τ)
end


function spectrum(ω, R::SymetricRTN)
    2 * R.b^2 * R.γ / (ω^2 + R.γ^2)
end


"""
$(TYPEDEF)

An ensemble of random telegraph noise.

$(FIELDS)
"""
struct EnsembleFluctuator{T}
    """A list of RTNs"""
    f::Vector{T}
end


function EnsembleFluctuator(b::AbstractArray{T}, ω::AbstractArray{T}) where T<:Number
    f = [SymetricRTN(x, y) for (x,y) in zip(b, ω)]
    EnsembleFluctuator(f)
end


function correlation(τ, E::EnsembleFluctuator)
    sum((x)->correlation(τ, x), E.f)
end


function spectrum(ω, E::EnsembleFluctuator)
    sum((x)->spectrum(ω, x), E.f)
end
