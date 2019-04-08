module QInterpolate

import Interpolations:interpolate,BSpline,Quadratic,Line,OnGrid,scale,gradient1,extrapolate
import Base.getindex

export construct_interpolations, COMPLEX_INTERPOLATION, COMPLEX_INTERPOLATION_ARRAY, REAL_INTERPOLATION_ARRAY, gradient, getindex

struct COMPLEX_INTERPOLATION
    real_part::AbstractArray{Float64,1}
    imag_part::AbstractArray{Float64,1}
end

struct COMPLEX_INTERPOLATION_ARRAY
    real_part::AbstractArray{AbstractArray{Float64, 1}, N} where N
    imag_part::AbstractArray{AbstractArray{Float64, 1}, N} where N
end

struct REAL_INTERPOLATION_ARRAY
    value::AbstractArray{AbstractArray{Float64, 1}, N} where N
end

function getindex(r_iter::REAL_INTERPOLATION_ARRAY, i...)
    getindex(r_iter.value, i...)
end

function getindex(c_iter::COMPLEX_INTERPOLATION_ARRAY, i...)
    r = getindex(c_iter.real_part, i...)
    i = getindex(c_iter.imag_part, i...)
    r,i
end

function (c_iter::COMPLEX_INTERPOLATION)(t)
    c_iter.real_part(t) + 1.0im*c_iter.imag_part(t)
end

function (c_iter::COMPLEX_INTERPOLATION_ARRAY)(t)
    res_size = size(c_iter.real_part)
    res = Array{ComplexF64}(undef, res_size)
    for i in eachindex(res)
        res[i] = c_iter.real_part[i](t) + 1.0im*c_iter.imag_part[i](t)
    end
    res
end

function (c_iter::REAL_INTERPOLATION_ARRAY)(t)
    res_size = size(c_iter.value)
    res = Array{Float64}(undef, res_size)
    for i in eachindex(res)
        res[i] = c_iter.value[i](t)
    end
    res
end

function construct_interpolations(x, y::Array{ComplexF64,1}; extrapolation="")
    yr = real(y)
    yi = imag(y)
    itpr = interpolate(yr, BSpline(Quadratic(Line(OnGrid()))))
    itpi = interpolate(yi, BSpline(Quadratic(Line(OnGrid()))))
    if lowercase(extrapolation) == "line"
        itpr = extrapolate(itpr, Line())
        itpi = extrapolate(itpi, Line())
    end
    sitpr = scale(itpr, x)
    sitpi = scale(itpi, x)
    COMPLEX_INTERPOLATION(sitpr, sitpi)
end

function construct_interpolations(x, y::Array{Float64,1}; extrapolation="")
    itp = interpolate(y,  BSpline(Quadratic(Line(OnGrid()))))
    if lowercase(extrapolation) == "line"
        itp = extrapolate(itp, Line())
    end
    sitp = scale(itp, x)
    return sitp
end

function construct_interpolations(x, y::Union{Array{Array{Float64,1},1},Array{Array{Float64,2},1}})
    res_size = size(y[1])
    interpolation_res = Array{AbstractArray}(undef,res_size)
    max_index = reduce(*,res_size)
    for j in range(1,length=max_index)
        interpolation_array = [y[i][j] for i in range(1,length=length(y))]
        itp = interpolate(interpolation_array,  BSpline(Quadratic(Line(OnGrid()))))
        interpolation_res[j] = scale(itp, x)
    end
    REAL_INTERPOLATION_ARRAY(interpolation_res)
end

function construct_interpolations(x, y::Union{Array{Array{ComplexF64,1},1},Array{Array{ComplexF64,2},1}})
    res_size = size(y[1])
    interpolation_res_real = Array{AbstractArray}(undef, res_size)
    interpolation_res_imag = Array{AbstractArray}(undef, res_size)
    max_index = reduce(*,res_size)
    for j in range(1, length=max_index)
        interpolation_array = [y[i][j] for i in range(1,length=length(y))]
        itpr = interpolate(real(interpolation_array), BSpline(Quadratic(Line(OnGrid()))))
        itpi = interpolate(imag(interpolation_array), BSpline(Quadratic(Line(OnGrid()))))
        interpolation_res_real[j] = scale(itpr, x)
        interpolation_res_imag[j] = scale(itpi, x)
    end
    COMPLEX_INTERPOLATION_ARRAY(interpolation_res_real, interpolation_res_imag)
end

function array_2_range(x)
    start = x[1]
    stop = x[end]
    range(start, stop=stop, length=length(x))
end

function gradient(c_iter::REAL_INTERPOLATION_ARRAY, t::Real)
    res = Array{Float64}(undef, size(c_iter.value))
    for i in eachindex(res)
        res[i] = gradient1(c_iter.value[i], t)
    end
    res
end

function gradient(x::COMPLEX_INTERPOLATION_ARRAY, t::Real)
    re = Array{ComplexF64}(undef, size(x.real_part))
    img = Array{ComplexF64}(undef, size(x.imag_part))
    for i in eachindex(re)
        re[i] = gradient1(x.real_part[i], t)
        img[i] = gradient1(x.imag_part[i], t)
    end
    re+1.0im*img
end

function gradient(x, t::Real)
    gradient1(x, t)
end

function gradient(x, t::AbstractArray{T, 1}) where T<:Real
    [gradient(x, i) for i in t]
end
end # module
