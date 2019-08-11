function prepare_u0(raw_u0, control = false)
    res = complex(raw_u0)
    if control == true
        if ndims(raw_u0) == 1
            res = DEVec(raw_u0, 1)
        else
            res = DEMat(raw_u0, 1)
        end
    end
    res
end


mutable struct DEVec{T} <: DEDataVector{T}
    x::Array{T,1}
    state::Int
end


mutable struct DEMat{T} <: DEDataMatrix{T}
    x::Array{T, 2}
    state::Int
end
