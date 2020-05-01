abstract type Manifold end

# fallback for out-of-place ops
retract(M::Manifold, x) = retract!(M, copy(x))
project_tangent(M::Manifold, g, x) = project_tangent!(M, copy(g), x)

# Fake objective function implementing a retraction
mutable struct ManifoldObjective{T<:NLSolversBase.AbstractObjective} <:
               NLSolversBase.AbstractObjective
    manifold::Manifold
    inner_obj::T
end
# TODO: is it safe here to call retract! and change x?
function NLSolversBase.value!(obj::ManifoldObjective, x)
    xin = retract(obj.manifold, x)
    value!(obj.inner_obj, xin)
end
function NLSolversBase.value(obj::ManifoldObjective)
    value(obj.inner_obj)
end
function NLSolversBase.gradient(obj::ManifoldObjective)
    gradient(obj.inner_obj)
end
function NLSolversBase.gradient(obj::ManifoldObjective, i::Int)
    gradient(obj.inner_obj, i)
end
function NLSolversBase.gradient!(obj::ManifoldObjective, x)
    xin = retract(obj.manifold, x)
    gradient!(obj.inner_obj, xin)
    project_tangent!(obj.manifold, gradient(obj.inner_obj), xin)
    return gradient(obj.inner_obj)
end
function NLSolversBase.value_gradient!(obj::ManifoldObjective, x)
    xin = retract(obj.manifold, x)
    value_gradient!(obj.inner_obj, xin)
    project_tangent!(obj.manifold, gradient(obj.inner_obj), xin)
    return value(obj.inner_obj)
end

"""Spherical manifold {||x|| = r}."""
struct Sphere{T} <: Manifold where {T<:Real}
    r::T
    Sphere(r::T) where {T<:Real} =
        r < 0 ? error("radius has to be a positive number!") : new{T}(r)
end
Sphere() = Sphere(1)
retract!(S::Sphere, x) = rmul!(x, S.r / norm(x))
project_tangent!(S::Sphere, g, x) = (g .-= (real(dot(x, g)) / S.r^2) .* x)


function ManifoldRetraction(M::Mfd; save = true) where {Mfd<:Manifold}
    affect!(integrator) = retract!(M, integrator.u)
    condition = (u, t, integrator) -> true
    save_positions = (false, save)
    DiscreteCallback(condition, affect!, save_positions = save_positions)
end
