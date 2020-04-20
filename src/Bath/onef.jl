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


correlation(τ, R::SymetricRTN) = R.b^2 * exp(-R.γ * τ)
spectrum(ω, R::SymetricRTN) = 2 * R.b^2 * R.γ / (ω^2 + R.γ^2)

construct_distribution(tf::Real, R::SymetricRTN) = Exponential(1 / R.γ / tf)
construct_distribution(tf::UnitTime, R::SymetricRTN) = Exponential(1 / R.γ)


"""
$(TYPEDEF)

An ensemble of random telegraph noise.

$(FIELDS)
"""
struct EnsembleFluctuator{T} <:AbstractBath
    """A list of RTNs"""
    f::Vector{T}
end


function EnsembleFluctuator(b::AbstractArray{T}, ω::AbstractArray{T}) where {T<:Number}
    f = [SymetricRTN(x, y) for (x, y) in zip(b, ω)]
    EnsembleFluctuator(f)
end

Base.length(E::EnsembleFluctuator) = Base.length(E.f)
correlation(τ, E::EnsembleFluctuator) = sum((x) -> correlation(τ, x), E.f)
spectrum(ω, E::EnsembleFluctuator) = sum((x) -> spectrum(ω, x), E.f)


construct_distribution(tf, E::EnsembleFluctuator) = product_distribution([construct_distribution(tf, x) for x in E.f])


function Base.show(io::IO, ::MIME"text/plain", E::EnsembleFluctuator)
    print(
        io,
        "Fluctuator ensemble with ",
        length(E),
        " fluctuators"
    )
end


function Base.show(io::IO, E::EnsembleFluctuator)
    print(
        io,
        "Fluctuator ensemble with ",
        length(E),
        " fluctuators"
    )
end
