"""
$(TYPEDEF)

Defines a fluctuator ensemble controller

# Fields

$(FIELDS)
"""
mutable struct FluctuatorControl{T} <: AbstractAnnealingControl
    """waitting time distribution for every fluctuators"""
    dist
    """cache for each fluctuator value"""
    b0
    """index of the fluctuator to be flipped next"""
    next_idx
    """time interval for next flip event"""
    next_τ
    """noise value"""
    n
end


function FluctuatorControl(tf, num::Int, E::EnsembleFluctuator)
    dist = construct_distribution(tf, E)
    b0 = [x.b for x in E.f] .* rand([-1, 1], length(dist), num)
    next_τ, next_idx = findmin(rand(dist, num))
    FluctuatorControl{num}(dist, b0, next_idx, next_τ, sum(b0, dims = 1)[:])
end

FluctuatorControl(tf, num, E, ::Nothing) = FluctuatorControl(tf, num, E)

(control::FluctuatorControl)() = control.n
get_controller_name(::FluctuatorControl) = :fluctuator_control


function reset!(ctrl::FluctuatorControl, ::DEFAULT_INITIALIZER_CLASS)
    ctrl.b0 =
        abs.(ctrl.b0) .* rand([-1, 1], length(ctrl.dist), size(ctrl.b0, 2))
    ctrl.n = sum(ctrl.b0, dims = 1)[:]
    next_state!(ctrl)
end


function reset!(ctrl::FluctuatorControl, initializer)
    ctrl.b0 =
        abs.(ctrl.b0) .* initializer((length(ctrl.dist), size(ctrl.b0, 2)))
    ctrl.n = sum(ctrl.b0, dims = 1)[:]
    next_state!(ctrl)
end

reset!(ctrl::FluctuatorControl, u0, initializer) = reset!(ctrl, initializer)


function build_callback(control::FluctuatorControl)
    fluctuator_time_choice = function (integrator)
        next_t = integrator.t + integrator.p.control.next_τ
        if next_t > integrator.sol.prob.tspan[2]
            return nothing
        else
            return next_t
        end
    end

    fluctuator_affect! = function (integrator)
        f = integrator.p.control
        f.n = sum(f.b0, dims = 1)[:]
        next_state!(f)
        u_modified!(integrator, false)
    end

    IterativeCallback(fluctuator_time_choice, fluctuator_affect!)
end


function build_callback(control::FluctuatorControl, controller_name::Symbol)
    fluctuator_time_choice = function (integrator)
        controller = getproperty(integrator.p.control, controller_name)
        next_t = integrator.t + controller.next_τ
        if next_t > integrator.sol.prob.tspan[2]
            return nothing
        else
            return next_t
        end
    end

    fluctuator_affect! = function (integrator)
        f = getproperty(integrator.p.control, controller_name)
        f.n = sum(f.b0, dims = 1)[:]
        next_state!(f)
        u_modified!(integrator, false)
    end

    IterativeCallback(fluctuator_time_choice, fluctuator_affect!)
end


function next_state!(f::FluctuatorControl)
    next_τ, next_idx = findmin(rand(f.dist, size(f.b0, 2)))
    f.next_τ = next_τ
    f.next_idx = next_idx
    f.b0[next_idx] *= -1
    nothing
end


"""
$(TYPEDEF)

Defines a fluctuator ensemble controller using DEDataArray

# Fields

$(FIELDS)
"""
mutable struct FluctuatorDEControl{T} <: AbstractAnnealingControl
    """Waitting time distribution for every fluctuators"""
    dist
    """Cache for each fluctuator value"""
    b0
    """Index of the fluctuator to be flipped next"""
    next_idx
    """Time interval for next flip event"""
    next_τ
    """Symbol used in DEDataArray"""
    sym::Symbol
end


function FluctuatorControl(tf, num::Int, E::EnsembleFluctuator, sym::Symbol)
    dist = construct_distribution(tf, E)
    b0 = [x.b for x in E.f] .* rand([-1, 1], length(dist), num)
    next_τ, next_idx = findmin(rand(dist, num))
    FluctuatorDEControl{num}(dist, b0, next_idx, next_τ, sym)
end


(f::FluctuatorDEControl)() = view(sum(f.b0, dims = 1), :)
get_controller_name(::FluctuatorDEControl) = :fluctuator_control


function reset!(ctrl::FluctuatorDEControl, u0, ::DEFAULT_INITIALIZER_CLASS)
    ctrl.b0 =
        abs.(ctrl.b0) .* rand([-1, 1], length(ctrl.dist), size(ctrl.b0, 2))
    setfield!(u0, ctrl.sym, Array(ctrl()))
    next_state!(ctrl)
end


function reset!(ctrl::FluctuatorDEControl, u0, initializer)
    ctrl.b0 =
        abs.(ctrl.b0) .* initializer((length(ctrl.dist), size(ctrl.b0, 2)))
    setfield!(u0, ctrl.sym, Array(ctrl()))
    next_state!(ctrl)
end


function build_callback(control::FluctuatorDEControl)
    fluctuator_time_choice = function (integrator)
        next_t = integrator.t + integrator.p.control.next_τ
        if next_t > integrator.sol.prob.tspan[2]
            return nothing
        else
            return next_t
        end
    end

    fluctuator_affect! = function (integrator)
        noise_value = Array(integrator.p.control())
        sym = getfield(integrator.p.control, :sym)
        for c in full_cache(integrator)
            setfield!(c, sym, noise_value)
#            c.n .= noise_value
        end
        next_state!(integrator.p.control)
        u_modified!(integrator, false)
    end

    IterativeCallback(fluctuator_time_choice, fluctuator_affect!)
end


function build_callback(control::FluctuatorDEControl, name::Symbol)
    fluctuator_time_choice = function (integrator)
        control = getproperty(integrator.p.control, name)
        next_t = integrator.t + control.next_τ
        if next_t > integrator.sol.prob.tspan[2]
            return nothing
        else
            return next_t
        end
    end

    fluctuator_affect! = function (integrator)
        control = getproperty(integrator.p.control, name)
        noise_value = Array(control())
        for c in full_cache(integrator)
            setfield!(c, control.sym, noise_value)
        end
        next_state!(control)
        u_modified!(integrator, false)
    end

    IterativeCallback(fluctuator_time_choice, fluctuator_affect!)
end


function next_state!(f::FluctuatorDEControl)
    next_τ, next_idx = findmin(rand(f.dist, size(f.b0, 2)))
    f.next_τ = next_τ
    f.next_idx = next_idx
    f.b0[next_idx] *= -1
    nothing
end


"""
$(TYPEDEF)

Defines stochastic system-bath coupling operator

# Fields

$(FIELDS)
"""
struct StochasticNoise <: AbstractOpenSys
    """System-bath coupling operator"""
    ops::AbstractCouplings
    """Symbol used for DEDataArray type"""
    sym
end


function (S::StochasticNoise)(A, n, tf::Real, t)
    A .+= -1.0im * tf * sum(n .* S.ops(t))
end


(S::StochasticNoise)(A, n, tf::UnitTime, t) = S(A, n, 1.0, t / tf)

function (S::StochasticNoise)(A, u::DEDataArray, tf::Real, t)
    A .+= -1.0im * tf * sum(getfield(u, S.sym) .* S.ops(t))
end


(S::StochasticNoise)(A, u::DEDataArray, tf::UnitTime, t) = S(A, u, 1.0, t / tf)
