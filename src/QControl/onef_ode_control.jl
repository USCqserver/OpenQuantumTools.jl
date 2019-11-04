"""
$(TYPEDEF)
DEDataVector type used for spin-fluctuator simulation.

# Fields
$(FIELDS)
"""
mutable struct DEFluctuatorVec{T} <: DEDataVector{T}
    """The state vector"""
    x::Array{T,1}
    """Noise value"""
    n::Float64
end


"""
$(TYPEDEF)
DEDataMatrix type used for spin-fluctuator simulation.

# Fields
$(FIELDS)
"""
mutable struct DEFluctuatorMat{T} <: DEDataMatrix{T}
    """The density matrix"""
    x::Array{T,2}
    """Noise value"""
    n::Float64
end
