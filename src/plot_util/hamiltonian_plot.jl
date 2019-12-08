@recipe function f(H::AbstractHamiltonian, s, lvl)
    y = []
    for x in s
        w ,_ = eigen_decomp(H, x; level=lvl)
        push!(y, w)
    end
    y = hcat(y...)'
    lab = ["E_$x" for x in (0:(lvl-1))']
    lab = [latexstring(x) for x in lab]
    label --> lab
    ylabel --> "GHz"
    xlabel --> "s"
    legend --> :best
    (s, y)
end
