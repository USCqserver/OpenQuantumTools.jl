@recipe function f(H::AbstractHamiltonian, s, lvl)
    y = []
    for x in s
        w ,_ = eigen_decomp(H, x; level=lvl)
        push!(y, w)
    end
    y = hcat(y...)' /2/π
    lab = ["E_$x" for x in 0:lvl]
    lab = [latexstring(x) for x in lab]
    label --> lab
    ylabel --> "GHz"
    xlabel --> "s"
    (s, y)
end
