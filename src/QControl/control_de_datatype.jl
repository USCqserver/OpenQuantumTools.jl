"""
$(TYPEDEF)

Base for types defining annealing controller that does not require DEDataArray in ODE solver.
"""
abstract type ParameterFreeControl <: AbstractAnnealingControl end


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

QTBase.check_positivity(x::Union{DEStateMachineMat, DENoiseMat}) = QTBase.check_positivity(x.x)
QTBase.check_unitary(x::Union{DEStateMachineMat, DENoiseMat}; rtol=1e-6, atol=1e-8) = QTBase.check_unitary(x.x, rtol=rtol, atol=atol)


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
sort_ctrs(ctrs, set::Nothing, args...) = sort_ctrs(ctrs, args...)

Base.length(c::ControlSet) = length(c.ctrs)
Base.iterate(c::ControlSet, state = 1) = Base.iterate(c.ctrs, state)


function build_controllers(bath)
end


function (h::DenseHamiltonian)(du, u::DEDataMatrix{T}, p::Real, t::Real) where T<:Complex
    fill!(du, 0.0+0.0im)
    H = h(t)
    LinearAlgebra.BLAS.gemm!('N', 'N', -1.0im * p, H, u.x, 1.0 + 0.0im, du.x)
    LinearAlgebra.BLAS.gemm!('N', 'N', 1.0im * p, u.x, H, 1.0 + 0.0im, du.x)
end
