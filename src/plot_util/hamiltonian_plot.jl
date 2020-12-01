@recipe function f(H::OpenQuantumBase.AbstractHamiltonian, s, lvl)
    y = []
    for x in s
        w, _ = eigen_decomp(H, x; lvl=lvl)
        push!(y, w)
    end
    y = hcat(y...)'
    lab = ["\$E_{$x}\$" for x in (0:lvl - 1)']
    label --> lab
    yguide --> "\$\\mathrm{GHz}\$"
    xguide --> "\$s\$"
    legend --> :best
    (s, y)
end
