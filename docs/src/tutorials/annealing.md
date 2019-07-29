# Annealing
This package implements [AbstractFactory](https://refactoring.guru/design-patterns/abstract-factory) pattern for potential quantum annealing process via an abstract type [`AbstractAnnealing`](@ref). A complete quantum annealing process is assembled from the following parts:
  1. Hamiltonian: Any object implements the `AbstractHamiltonian` interface
  2. Initial state: A state vector/density matrix
  3. (Optional)System bath coupling -- system part
  4. (Optional)System bath coupling -- bath part
  5. (Optional)Additional control protocols
