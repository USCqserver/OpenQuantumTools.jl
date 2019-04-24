"""
    check_unitary(𝐔; rtol=1e-6, atol=1e-8)

Test if `𝐔` is a unitary matrix. The function checks how close both `` 𝐔𝐔^† `` and `` 𝐔^†𝐔 `` are to `` I ``, with relative and absolute error given by `rtol`, `atol`.

# Examples
```julia-repl
julia> check_unitary(exp(-1.0im*5*0.5*σx))
true
```
"""
function check_unitary(𝐔::Array{T,2}; rtol=1e-6, atol=1e-8) where T<:Number
    a1 = isapprox(𝐔*𝐔', Matrix{eltype(𝐔)}(I,size(𝐔)),rtol=rtol,atol=atol)
    a2 = isapprox(𝐔'*𝐔, Matrix{eltype(𝐔)}(I,size(𝐔)),rtol=rtol,atol=atol)
    a1 && a2
end
