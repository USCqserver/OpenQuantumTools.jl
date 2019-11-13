"""
$(TYPEDEF)
DEDataVector type for finite state machine control

# Fields
$(FIELDS)
"""
mutable struct DEStateMachineVec{T} <: DEDataVector{T}
    """DiffEq vector"""
    x::Array{T,1}
    """Current state"""
    state::Int
end


"""
$(TYPEDEF)
DEDataMatrix type for finite state machine control

# Fields
$(FIELDS)
"""
mutable struct DEStateMachineMat{T} <: DEDataMatrix{T}
    """DiffEq matrix"""
    x::Array{T,2}
    """Current state"""
    state::Int
end


"""
$(TYPEDEF)
DEDataVector type for noise injection

# Fields
$(FIELDS)
"""
mutable struct DENoiseVec{T} <: DEDataVector{T}
    """DiffEq vector"""
    x::Array{T,1}
    """Current noise value"""
    n::Float64
end


"""
$(TYPEDEF)
DEDataMatrix type for noise injection

# Fields
$(FIELDS)
"""
mutable struct DENoiseMat{T} <: DEDataMatrix{T}
    """DiffEq matrix"""
    x::Array{T,2}
    """Current noise value"""
    n::Float64
end


adjust_u0_with_control(u0, ::Nothing) = u0
construct_callback(::Nothing, ::Symbol) = nothing


function (h::DenseHamiltonian)(du, u::DEDataMatrix{T}, p::Real, t::Real) where T<:Complex
    fill!(du, 0.0+0.0im)
    H = h(t)
    LinearAlgebra.BLAS.gemm!('N', 'N', -1.0im * p, H, u.x, 1.0 + 0.0im, du.x)
    LinearAlgebra.BLAS.gemm!('N', 'N', 1.0im * p, u.x, H, 1.0 + 0.0im, du.x)
end
