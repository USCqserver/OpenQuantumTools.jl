# Bath Module
## Ohmic Bath Object
```@docs
OhmicBath
Ohmic(η, fc, T)
correlation(τ, params::OhmicBath)
polaron_correlation(τ, params::OhmicBath)
γ(w, params::OhmicBath)
S(w, params::OhmicBath; atol=1e-7)
interpolate_spectral_density(ω_grid::AbstractRange{T}, params::OhmicBath) where T<:Number
```
## Hybrid Ohmic Bath
```@docs
HybridOhmicBath
HybridOhmic(W, η, fc, T)
reorganization_energy(bath::HybridOhmicBath)
ohmic_correlation(τ, bath::HybridOhmicBath[, a=1])
polaron_correlation(τ, bath::HybridOhmicBath[, a=1])
GH(ω, bath::HybridOhmicBath[, a=1])
```
