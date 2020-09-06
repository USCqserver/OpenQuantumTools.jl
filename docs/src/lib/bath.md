# Bath Module
## Ohmic bath
```@docs
OhmicBath
Ohmic
```
## Common interface
```@docs
correlation(τ, bath::OhmicBath)
γ(ω, bath::OhmicBath)
spectrum(ω, bath::OhmicBath)
S(w, bath::OhmicBath; atol=1e-7)
```