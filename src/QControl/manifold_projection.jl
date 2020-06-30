abstract type Manifold end

# fallback for out-of-place ops
retract(M::Manifold, x) = retract!(M, copy(x))
project_tangent(M::Manifold, g, x) = project_tangent!(M, copy(g), x)

"""Spherical manifold {||x|| = r}."""
struct Sphere{T} <: Manifold where {T<:Real}
    r::T
    Sphere(r::T) where {T<:Real} =
        r < 0 ? error("radius has to be a positive number!") : new{T}(r)
end
Sphere() = Sphere(1)
retract!(S::Sphere, x) = rmul!(x, S.r / norm(x))
#project_tangent!(S::Sphere, g, x) = (g .-= (real(dot(x, g)) / S.r^2) .* x)

function ManifoldRetraction(M::Mfd; save = true) where {Mfd<:Manifold}
    affect!(integrator) = retract!(M, integrator.u)
    condition = (u, t, integrator) -> true
    DiscreteCallback(condition, affect!, save_positions = (false, save))
end
