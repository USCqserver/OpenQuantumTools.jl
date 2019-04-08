module Integration

using Reexport
@reexport using QuadGK
#@reexport using HCubature
#@reexport using Cuba

export integrate_1d

include("scipy_integration.jl")
#include("cuba_integration.jl")
#include("hcubature_integration.jl")
include("quadgk_integration.jl")

"""
    integrate_1d(f, a, b; args=(), rtol=1e-6, atol=1e-8)

Calculate numerical integration of function `f` from `a` to `b`.
"""
function integrate_1d(f, a, b; args=(), rtol=1e-6, atol=1e-8, kwargs...)
    integrand, lower_limit, upper_limit = integral_limit_transfromation(f, a, b, args)
    quadgk(integrand, lower_limit, upper_limit, rtol=rtol, atol=atol; kwargs...)
end

end # module
