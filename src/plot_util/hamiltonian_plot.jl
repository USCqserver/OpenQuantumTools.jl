@recipe function f(H::AbstractHamiltonian, s, lvl; tol=1e-8)
    y = []
    for x in s
        w ,_ = eigen_decomp(H, x; level=lvl, tol=tol)
        push!(y, w)
    end
    y = hcat(y...)'
    lab = ["E_$x" for x in (0:(lvl-1))']
    lab = [latexstring(x) for x in lab]
    label --> lab
    yguide --> L"\mathrm{GHz}"
    xguide --> L"s"
    legend --> :best
    (s, y)
end
