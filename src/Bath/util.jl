build_correlation(bath::AbstractBath, tf::Real) =
    (s) -> correlation(s * tf, bath)
build_correlation(bath::AbstractBath, ::UnitTime) = (τ) -> correlation(τ, bath)
build_spectrum(bath::AbstractBath) = (ω) -> spectrum(ω, bath)

function build_redfield(
    coupling::AbstractCouplings,
    unitary,
    tf::Union{Real,UnitTime},
    bath::AbstractBath;
    atol = 1e-8,
    rtol = 1e-6,
)
    cfun = build_correlation(bath, tf)
    Redfield(coupling, unitary, cfun, atol = atol, rtol = rtol)
end

function build_davies(
    coupling::AbstractCouplings,
    bath::AbstractBath,
    ω_range,
    lambshift::Bool,
)
    if lambshift == true
        if isempty(ω_range)
            S_loc = (ω) -> S(ω, bath)
        else
            s_list = [S(ω, bath) for ω in ω_range]
            S_loc = construct_interpolations(ω_range, s_list)
        end
    else
        S_loc = (ω) -> 0.0
    end
    DaviesGenerator(coupling, build_spectrum(bath), S_loc)
end

function build_CGME(
    coupling::AbstractCouplings,
    unitary,
    tf::Union{Real,UnitTime},
    bath::AbstractBath;
    atol = 1e-8,
    rtol = 1e-6,
    Ta = nothing,
)
    Ta = Ta == nothing ? coarse_grain_timescale(bath, tf) : Ta
    cfun = build_correlation(bath, tf)
    CGOP(coupling, unitary, cfun, Ta, atol = atol, rtol = rtol)
end

"""
    lambshift(w, γ; atol=1e-7)

Calculate the Lamb shift of spectrum `γ`. `atol` is the absolute tolerance for Cauchy principal value integral.
"""
function lambshift(w, γ; atol = 1e-7)
    g(x) = γ(x) / (x - w)
    cpv, cperr = cpvagk(γ, w, w - 1.0, w + 1.0)
    negv, negerr = quadgk(g, -Inf, w - 1.0)
    posv, poserr = quadgk(g, w + 1.0, Inf)
    v = cpv + negv + posv
    err = cperr + negerr + poserr
    if (err > atol) || (isnan(err))
        @warn "Absolute error of integration is larger than the tolerance."
    end
    -v / 2 / pi
end

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

function coarse_grain_timescale(
    bath::AbstractBath,
    lim;
    rtol = sqrt(eps()),
    atol = 0,
)
    τsb, err_sb = τ_SB(bath, rtol = rtol, atol = atol)
    τb, err_b = τ_B(bath, lim, τsb, rtol = rtol, atol = atol)
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
