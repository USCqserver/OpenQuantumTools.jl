@recipe function f(H::AbstractHamiltonian, s, lvl, eig_init = EIGEN_DEFAULT)
    y = []
    for x in s
        w, _ = eigen_decomp(H, x; lvl = lvl, eig_init = EIGEN_DEFAULT)
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
