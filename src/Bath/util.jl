"""
    function τ_SB(cfun; lim=Inf, rtol=sqrt(eps()), atol=0)

Calculate τ_SB of the bath correlation function. It is defined as the integration of the absolute value of bath correlation function from zero to infinity. `lim` is the upper limit of the integration. `atol` and `rtol` are the absolute and relative error of the integration.
"""
function τ_SB(cfun; lim=Inf, rtol=sqrt(eps()), atol=0)
    res, err = quadgk((x)->abs(cfun(x)), 0, lim, rtol=rtol, atol=atol)
    1/res, err
end


"""
    function τ_B(cfun, lim, τsb; rtol=sqrt(eps()), atol=0)

Calculate the bath correlation time τ_B. The upper limit `lim` and `τsb` need to be manually specified. `atol` and `rtol` are the absolute and relative error for the integration.
"""
function τ_B(cfun, lim, τsb; rtol=sqrt(eps()), atol=0)
    res, err = quadgk((x)->x*abs(cfun(x)), 0, lim, rtol=rtol, atol=atol)
    res*τsb, err*τsb
end
