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
