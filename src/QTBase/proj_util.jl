function _init_proj(H::SparseMatrixCSC{T, V}, dH::SparseMatrixCSC{T, V}, lvl, tol) where T<:Number where V<:Int
    w, v = eigs(H, nev=lvl, which=:SR, tol=tol)
end

function _proj_step!(ref, H::SparseMatrixCSC{T, V}, dH::SparseMatrixCSC{T, V}, lvl, tol) where T<:Number where V<:Int
    # note this function also return a vector w, this is not a good practive but I need to think more for an elegant solution
    w, v = eigs(H, nev=lvl, which=:SR, tol=tol, v0=ref[:,1])
    for i in 1:lvl
        if view(v, :, i)' * view(ref, :, i) < 0
            ref[:, i] = -v[:, i]
        else
            ref[:, i] = v[:, i]
        end
    end
    w
end

function _calc_angle(dH, v)
end

function proj_low_lvl(hfun, dhfun, interaction, s_axis::AbstractArray{T, 1}; ref=nothing, tol=1e-4, lvl=2) where T<:Number
    ev = Array{Array{Float64, 1}, 1}(undef, length(s_axis))
    dθ = Array{Array{Float64, 2}, 1}(undef, length(s_axis))
    op = Array{Array{Array{Float64, 2},1},1}(undef, length(s_axis))

    H = hfun(s_axis[1])
    dH = dhfun(s_axis[1])
    if ref==nothing
        w, v = _init_proj(H, dH, lvl, tol)
    else
        w = _proj_step(ref, H, dH, lvl, tol)
    end
    push!(ev, w)

end
