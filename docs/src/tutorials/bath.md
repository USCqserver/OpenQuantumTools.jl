# Bath
We start with a general form of system bath interaction Hamiltonian
```math
  H_{SB} = \sum_\alpha S_\alpha\otimes B_\alpha
```
and define each $B_\alpha$ (together with the bath Hamiltonian $H_B$) as an `AbstractBath` object. Currently, there are two build-in concrete types of this class -- [`OhmicBath`](@ref) and [`HybridOhmicBath`](@ref).
## Ohmic Bath
## Hybrid Ohmic Bath
A good reference paper for the hybrid ohmic bath model is [Lanting et al.](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.83.180502) The model is generalized in [Smirnov and Amin](https://iopscience.iop.org/article/10.1088/1367-2630/aae79c/meta). The basic idea is, the noise spectrum of this bath model can be split into the low frequency and high frequency parts
```math
  S(\omega) = S^L(ω) + S^H(ω)
```
The high frequency part of the bath is the Ohmic and the low frequency part can be something like the $1/f$ noise. The benefit of this formalism is that, instead of the entire spectrum density, $S^L(\omega)$ can be parametrized by a single parameter in the macroscopic resonant tunneling(MRT) experiment. 
