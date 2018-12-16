using DifferentialEquations
export calculate_unitary, unitary_check

function calculate_unitary(𝐇;rtol=1e-6,atol=1e-8)
    u0 = Matrix{ComplexF64}(I, size(𝐇(0.0)))
    function f(du, u, p, t)
        hmat = -1.0im * 𝐇(t)
        mul!(du, hmat, u)
    end
    prob = ODEProblem(f, u0, (0.0,1.0))
    sol = solve(prob,Tsit5(),reltol=rtol,abstol=atol)
end

"""
    unitary_check(𝐔; rtol=1e-6, atol=1e-8)

Test if `𝐔` is a unitary matrix. The function checks how close both ``` 𝐔*𝐔' ``` and ``` 𝐔'*𝐔 ``` are to ``` I ```, with relative and absolute error given by `rtol`, `atol`.

# Examples
```julia-repl
julia> unitary_check(exp(-1.0im*5*0.5*σx))
true
```
"""
function unitary_check(𝐔::Array{T,2}; rtol=1e-6, atol=1e-8) where T<:Number
    a1 = isapprox(𝐔*𝐔', Matrix{eltype(𝐔)}(I,size(𝐔)),rtol=rtol,atol=atol)
    a2 = isapprox(𝐔'*𝐔, Matrix{eltype(𝐔)}(I,size(𝐔)),rtol=rtol,atol=atol)
    a1 && a2
end
