struct LinearOperator
    " List of time dependent functions "
    f::Array{Function, 1}
    " List of constant matrices "
    m::Array{Array{T,2},1} where T<:Number
end

function (h::LinearOperator)(t::Real)
    res = zeros(eltype(h.m[1]), size(h.m[1]))
    for (f,m) in zip(h.f,h.m)
        axpy!(f(t),m,res)
    end
    res
end

function multiply(a, h::LinearOperator)
    LinearOperator(h.f, a*h.m)
end

function multiply!(h::LinearOperator, a)
    for i in eachindex(h.m)
        h.m[i] *= a
    end
end

function update!(u, t::Real, h::LinearOperator)
    for (f,m) in zip(h.f,h.m)
        axpy!(f(t),m,u)
    end
end

function comm!(du, u, t::Real, h::LinearOperator; p = -1.0im)
    for (f, m) in zip(h.f, h.m)
        c = p * f(t)
        gemm!('N', 'N', c, m, u, 1.0, du)
        gemm!('N', 'N', -c, u, m, 1.0, du)
    end
end
