export monte_carlo_integrate, monte_carlo_𝒯integrate

function monte_carlo_integrate(func, a::Array{Float64,1}, b::Array{Float64,1}, ncomp::Int; args=(), method::Symbol=:vegas, kwargs...)
    integrand = generate_integrand(func, a, b, args)
    ndim = length(a)
    @eval begin
        $(method)($integrand,$ndim,$ncomp; $kwargs...)
    end
end

function generate_integrand(func, a, b, args)
    (x,f)-> f .= func(a + (b-a).*x, args...)
end

function monte_carlo_𝒯integrate(func, t, ndim::Int, ncomp::Int; args=(), method::Symbol=:vegas)
    integrand = generate_𝒯integrand_cuba(func, t, args=args)
    @eval begin
        $(method)($integrand,$ndim,$ncomp)
    end
end

function generate_𝒯integrand_cuba(func, t; args=())
    function integrand(x, f)
        y = 1.0
        x[1] = x[1]*t
        for i in 2:length(x)
            x[i] = x[i]*x[i-1]
            y = y*x[i-1]
        end
        f.=func(x, args...)*y*t
    end
end
