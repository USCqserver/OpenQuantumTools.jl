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


function Spectrum(ω, R::SymetricRTN)
    2 * R.b^2 * R.γ / (ω^2 + R.γ^2)
end
