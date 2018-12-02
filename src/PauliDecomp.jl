export interaction_decompose

function interaction_decompose(unitary, operator, coordinate)
    dim = size(coordinate[1])[1]
    t_dim = length(unitary)
    res = Array{Array{Float64, 1},1}(undef, t_dim)
    for i in eachindex(unitary)
        u_op = unitary[i]' * operator * unitary[i]
        c_vec = [tr(u_op*c)/dim for c in coordinate]
        res[i] = real(c_vec)
    end
    res
end
