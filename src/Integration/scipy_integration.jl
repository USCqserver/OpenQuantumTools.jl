using PyCall

export scipy_quad, cauchy_principal_value

const scipy_integrate = PyNULL()

function __init__()
    copy!(scipy_integrate, pyimport_conda("scipy.integrate", "scipy"))
end

function scipy_quad(func, a, b; args=(), full_output=0, atol=1.49e-8, rtol=1.49e-8, limit=50, points=nothing, weight=nothing, wvar=nothing, wopts=nothing, maxp1=50, limlst=50)
    scipy_integrate[:quad](func, a, b; args=args, full_output=full_output, epsabs=atol, epsrel=rtol, limit=limit, points=points, weight=weight, wvar=wvar, wopts=wopts, maxp1=maxp1, limlst=limlst)
end

function cauchy_principal_value(func, pole; atol = 1e-8)
    s = scipy_integrate[:quad](func,pole-1.0,pole+1.0,weight="cauchy",wvar=pole)
    g(x) = func(x) / (x - pole)
    n = scipy_integrate[:quad](g,-Inf,pole-1.0)
    p = scipy_integrate[:quad](g,pole+1.0,Inf)
    err = n[2]+s[2]+p[2]
    if (err > atol) || (isnan(err))
        @warn "Absolute error of integration is larger than the tolerance."
    end
    s[1]+n[1]+p[1] , err
end
