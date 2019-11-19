"""
$(TYPEDEF)
DEDataArray type for finite state machine control.

# Fields
$(FIELDS)
"""
mutable struct DEStateMachineArray{T, N} <: DEDataArray{T, N}
    """Array data"""
    x::Array{T, N}
    """Current state"""
    state::Int
end


"""
$(TYPEDEF)
DEDataArray type for noise injection.

# Fields
$(FIELDS)
"""
mutable struct DENoiseArray{T, N} <: DEDataArray{T, N}
    """Array data"""
    x::Array{T, N}
    """Current noise value"""
    n::Vector{Float64}
end


const DEStateMachineVec{T} = DEStateMachineArray{T, 1}
const DEStateMachineMat{T} = DEStateMachineArray{T, 2}
const DENoiseVec{T} = DENoiseArray{T, 1}
const DENoiseMat{T} = DENoiseArray{T, 2}


"""
$(TYPEDEF)
DEDataArray type for both finite state machine control and noise injection.

# Fields
$(FIELDS)
"""
mutable struct DESTNoiseArray{T, N} <: DEDataArray{T, N}
    """Array data"""
    x::Array{T, N}
    """Current state"""
    state::Int
    """Current noise value"""
    n::Vector{Float64}
end


const DESTNoiseVec{T} = DESTNoiseArray{T, 1}
const DESTNoiseMat{T} = DESTNoiseArray{T, 2}


"""
$(TYPEDEF)
ControlSet to combine multiple control protocols.

# Fields
$(FIELDS)
"""
struct ControlSet{T<:Tuple} <: AbstractAnnealingControl
    """Tuple for control objects"""
    ctrs::T
end


ControlSet(::Nothing) = nothing
ControlSet(ctrs::Union{AbstractAnnealingControl, Nothing}...) = ControlSet(sort_ctrs((),ctrs...))

sort_ctrs(ctrs) = ctrs
sort_ctrs(ctrs, ctr::AbstractAnnealingControl, args...) = sort_ctrs((ctrs..., ctr), args...)
sort_ctrs(ctrs, set::ControlSet, args...) = sort_ctrs((ctrs..., set.ctrs...), args...)


function (h::DenseHamiltonian)(du, u::DEDataMatrix{T}, p::Real, t::Real) where T<:Complex
    fill!(du, 0.0+0.0im)
    H = h(t)
    LinearAlgebra.BLAS.gemm!('N', 'N', -1.0im * p, H, u.x, 1.0 + 0.0im, du.x)
    LinearAlgebra.BLAS.gemm!('N', 'N', 1.0im * p, u.x, H, 1.0 + 0.0im, du.x)
end
