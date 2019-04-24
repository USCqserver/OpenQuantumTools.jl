using QTool, Test

function gamma(x)
    if x==0
        return 1
    elseif x>0
        return x+1
    else
        return (1-x)*exp(x)
    end
end

function Sf(x)
    return x+0.1
end

H = -0.5 * standard_driver(2) + 0.5*(0.1*σz⊗σi + 0.5*σz⊗σz)
op = σz⊗σi + σi⊗σz
w, v = eigen(Hermitian(H))
state = (v[:, 1] + v[:, 2] + v[:, 3])/sqrt(3)
rho = state* state'
u = v'*rho*v
L_ij = Array{Array{Complex{Float64},2}, 2}(undef, (4, 4))
A_ij = Array{Complex{Float64}, 2}(undef, (4, 4))
for i in 1:4
    for j in 1:4
        A_ij[i,j] = v[:, i]' * op * v[:, j]
        L_ij[i,j] = v[:, i]' * op * v[:, j] * v[:, i] * v[:, j]'
    end
end
drho = []
hls= []
w_ab = Array{Float64, 2}(undef, (4,4))
for i in 1:4
    for j in 1:4
        w_ab[i,j] = w[j] - w[i]
        if i!=j
            T = L_ij[i,j]' * L_ij[i,j]
            push!(drho, gamma(w[j]-w[i])*(L_ij[i,j]*rho*L_ij[i,j]'-0.5*(T*rho+rho*T)))
            push!(hls, T*Sf(w[j]-w[i]))
        end
        T = L_ij[i,i]' * L_ij[j,j]
        push!(drho, gamma(0)*(L_ij[i,i]*rho*L_ij[j,j]'-0.5*(T*rho+rho*T)))
        push!(hls, T*Sf(0))
    end
end
drho = sum(drho)
hls = sum(hls)
drho = drho -1.0im * (hls*rho - rho*hls)
gm = gamma.(w_ab)
sm = Sf.(w_ab)
du = zeros(ComplexF64, (4,4))
adiabatic_me_update!(du, u, A_ij, gm, sm)
@test isapprox(v * du * v', drho, atol=1e-6, rtol=1e-6)
