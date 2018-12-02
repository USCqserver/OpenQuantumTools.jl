export eigen_value_eval, eigen_state_eval, inst_population

"""
    eigen_state_eval(hfun, t[, levels=[1,]])

Calculate the eigen states of Hamiltonian 'hfun' at each points of t. The output levels are specified by 'levels' argument.
"""
function eigen_state_eval(hfun, t::AbstractArray{Float64,1}; levels::Array{Int64,1}=[1,])
    res_dim = length(levels)
    res = Array{Array{Array{ComplexF64,1},1},1}(undef, length(t))
    for i in eachindex(t)
        eig_sys = eigen(Hermitian(hfun(t[i])))
        res[i] = [eig_sys.vectors[:,x] for x in levels]
    end
    for i in range(2,stop=length(t))
        for j in range(1,length=res_dim)
            if real(res[i][j]'*res[i-1][j]) < 0
                res[i][j] = -res[i][j]
            end
        end
    end
    res
end

"""
    eigen_value_eval(hfun, t[, levels])

Calculate the eigen values of Hamiltonian 'hfun' at each points of 't'. The output levels are specified by 'levels' argument.
"""
function eigen_value_eval(hfun, t::AbstractArray{Float64,1}; levels::Array{Int64,1}=[1,])
    res_dim = length(levels)
    res = Array{Array{Float64,1},1}(undef, length(t))
    for i in eachindex(t)
        eig_sys = eigen(Hermitian(hfun(t[i])))
        res[i] = [eig_sys.values[x] for x in levels]
    end
    res
end

"""
    function inst_population(ode_sol, hamiltonian; level::Int64=1)

Calculate the population of instantaneous eigen states specified by 'level', given array of quantum states 'ode_sol'. 'ode_sol' should have fields: 't' and 'u'.
"""
function inst_population(ode_sol, hamiltonian; level::Int64=1)
    pop = Array{Array{Float64, 1}, 1}(undef, length(ode_sol.t))
    for i in eachindex(ode_sol.t)
        hmat = hamiltonian(ode_sol.t[i])
        eig_sys = eigen(Hermitian(hmat))
        if ndims(ode_sol.u[i]) == 1
            inst_state = view(eig_sys.vectors,:,1:level)'
            pop[i] = abs2.(inst_state * ode_sol.u[i])
        elseif ndims(ode_sol.u[i]) == 2
            temp = Array{Float64, 1}(undef, level)
            for j in range(1, length=level)
                inst_state = view(eig_sys.vectors,:,j)
                temp[j] = real(inst_state'*ode_sol.u[i]*inst_state)
            end
            pop[i] = temp
        end
    end
    if level==1
        return [x[1] for x in pop]
    else
        return pop
    end
end
