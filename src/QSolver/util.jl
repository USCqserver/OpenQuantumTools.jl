"""
$(SIGNATURES)

Prepare the initial state in proper type and dimension for ODE solvers. `type` specifies the dimension of the initial state: `:v` is 1-D state vector and `:m` is 2-D density matrix. `vectorize` indicate whether to vectorize the density matrix.
"""
function build_u0(raw_u0, type; vectorize::Bool = false)
    res = complex(raw_u0)
    if type == :v && ndims(res) != 1
        throw(ArgumentError("Cannot convert density matrix to state vector."))
    elseif type == :m && ndims(res) == 1
        res = res * res'
    elseif ndims(res) < 1 || ndims(res) > 2
        throw(ArgumentError("u0 can either be a vector or matrix."))
    end
    if vectorize == true
        res = res[:]
    end
    # use StaticArrays for initial state will cause some method to fail
    # if (ndims(res) == 1) && (length(res) <= 10)
    #     res = OpenQuantumBase.MVector{length(res)}(res)
    # elseif (ndims(res) == 2) && (size(res, 1) <= 10)
    #     res = OpenQuantumBase.MMatrix{size(res, 1),size(res, 2)}(res)
    # end
    res
end

function vectorize_cache(cache)
    # use StaticArray if the dimension is less than 3
    if size(cache, 1) <= 3
        OpenQuantumBase.@MMatrix zeros(eltype(cache), size(cache, 1)^2, size(cache, 2)^2)
    else
        one(cache) ⊗ cache
    end
end
