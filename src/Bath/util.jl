"""
    function τ_SB(cfun; lim=Inf, rtol=sqrt(eps()), atol=0)

Calculate τ_SB of the bath correlation function. It is defined as the integration of the absolute value of bath correlation function from zero to infinity. `lim` is the upper limit of the integration. `atol` and `rtol` are the absolute and relative error of the integration.
"""
function τ_SB(cfun; lim = Inf, rtol = sqrt(eps()), atol = 0)
    res, err = quadgk((x) -> abs(cfun(x)), 0, lim, rtol = rtol, atol = atol)
    1 / res, err / abs2(res)
end


τ_SB(bath::AbstractBath; lim = Inf, rtol = sqrt(eps()), atol = 0) =
    τ_SB((t) -> correlation(t, bath), lim = lim, rtol = rtol, atol = atol)


"""
    function τ_B(cfun, lim, τsb; rtol=sqrt(eps()), atol=0)

Calculate the bath correlation time τ_B. The upper limit `lim` and `τsb` need to be manually specified. `atol` and `rtol` are the absolute and relative error for the integration.
"""
function τ_B(cfun, lim, τsb; rtol = sqrt(eps()), atol = 0)
    res, err = quadgk((x) -> x * abs(cfun(x)), 0, lim, rtol = rtol, atol = atol)
    res * τsb, err * τsb
end


τ_B(bath::AbstractBath, lim, τsb; rtol = sqrt(eps()), atol = 0) =
    τ_B((t) -> correlation(t, bath), lim, τsb; rtol = rtol, atol = atol)


function coarse_grain_timescale(bath::AbstractBath, lim; rtol=sqrt(eps()), atol=0)
    τsb, err_sb = τ_SB(bath, rtol = rtol, atol = atol)
    τb, err_b = τ_B(bath, lim, τsb, rtol=rtol, atol = atol)
    sqrt(τsb * τb / 5)
end


function g_ULE(t, spectrum; rtol = sqrt(eps()), atol = 0)
    res, err = quadgk(
        (w) -> sqrt(spectrum(w)) * exp(-1.0im * t * w) / sqrt(2π),
        -Inf,
        Inf,
        rtol = rtol,
        atol = atol,
    )
end

g_ULE(t, bath::AbstractBath; rtol = sqrt(eps()), atol = 0) =
    g_ULE(t, (w) -> spectrum(w, bath), rtol = rtol, atol = atol)
