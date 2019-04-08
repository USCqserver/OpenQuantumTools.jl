function semi_infinite_transform(f, a, b, args)
    if b == Inf
        return (y)->f(a+y/(1-y), args...)/(1-y)^2, 0.0, 1.0
    elseif b == -Inf
        return (y)->f(a+(1-y)/y, args...)/y^2, 1.0, 0.0
    elseif a == Inf
        return (y)->f(b+y/(1-y), args...)/(1-y)^2, 1.0, 0.0
    elseif a == -Inf
        return (y)->f(b+(1-y)/y, args...)/y^2, 0.0, 1.0
    else
        @error "Invalid integration limit"
    end
end

function infinite_transformation(f, a, b, args)
    if a == b
        @error "Invalid integration limit"
    else
        lower_limit = a < 0 ? 0.0 : 1.0
        upper_limit = 1.0 - lower_limit
    end
    (y)->f((2*y-1)/(1-y)/y, args...)*(2*y^2-2*y+1)/(1-y)^2/y^2, lower_limit, upper_limit
end

function integral_limit_transfromation(f, a, b, args)
    if isfinite(a) && isfinite(b)
        return (y)->f(y, args...), a, b
    elseif !isfinite(a) && !isfinite(b)
        return infinite_transformation(f, a, b, args)
    else
        return semi_infinite_transform(f, a, b, args)
    end
end
