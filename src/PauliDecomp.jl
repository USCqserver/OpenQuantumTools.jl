export matrix_decompose

"""
    matrix_decompose(mat::Array{T,2}, basis::Array{Array{T,2},1})

Decompse matrix `mat` onto matrix basis `basis`
"""
function matrix_decompose(mat::Array{T,2}, basis::Array{Array{T,2},1}) where T<:Number
    dim = size(basis[1])[1]
    [tr(mat*b)/dim for b in basis]
end
