# Bath Module
## Ohmic bath
```@docs
OhmicBath
Ohmic
```
## Hybrid Ohmic bath
```@docs
HybridOhmicBath
HybridOhmic
```
## Spin-fluctuator
```@docs
EnsembleFluctuator
```
## Others
```@docs
CustomBath
CorrelatedBath
```
## Common interface
```@docs
correlation(τ, bath::OhmicBath)
γ(ω, bath::OhmicBath)
spectrum(ω, bath::OhmicBath)
S(w, bath::OhmicBath; atol=1e-7)
```