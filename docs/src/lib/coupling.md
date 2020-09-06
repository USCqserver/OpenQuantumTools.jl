# Couplings
## Constant couplings
```@docs
ConstantCouplings
ConstantCouplings(mats::Union{Vector{Matrix{T}},Vector{SparseMatrixCSC{T,Int}}}; unit=:h) where {T<:Number}
ConstantCouplings(c::Vector{T}; sp = false, unit=:h) where T <: AbstractString
collective_coupling(op, num_qubit; sp = false, unit = :h)
```
## Time dependent couplings
```@docs
TimeDependentCoupling
TimeDependentCouplings
```