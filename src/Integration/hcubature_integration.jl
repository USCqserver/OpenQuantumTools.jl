export 𝒯hcubature, generate_𝒯integrand

function generate_𝒯integrand(func, t; args=())
    function integrand(x)
        y = 1.0
        sv = similar(x)
        sv[1] = x[1]*t
        for i in 2:length(x)
            sv[i] = x[i]*sv[i-1]
            y = y*sv[i-1]
        end
        func(sv, args...)*y*t
    end
end

function 𝒯hcubature(f, t, ndim, args=(); rtol=1e-6, atol=1e-8)
    integrand = generate_𝒯integrand(f, t, args=args)
    a = zeros(ndim)
    b = ones(ndim)
    hcubature(integrand, a, b, rtol=rtol, atol=atol)
end
