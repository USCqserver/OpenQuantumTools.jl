"""
$(TYPEDEF)

Defines a fluctuator ensemble controller

# Fields

$(FIELDS)
"""
mutable struct FluctuatorControl{T} <: AbstractAnnealingControl
    """waitting time distribution for every fluctuators"""
    dist::Any
    """cache for each fluctuator value"""
    b0::Any
    """index of the fluctuator to be flipped next"""
    next_idx::Any
    """time interval for next flip event"""
    next_τ::Any
    """noise value"""
    n::Any
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

    IterativeCallback(
        fluctuator_time_choice,
        fluctuator_affect!,
        save_positions = (false, true),
    )
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

    IterativeCallback(
        fluctuator_time_choice,
        fluctuator_affect!,
        save_positions = (false, true),
    )
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
    dist::Any
    """Cache for each fluctuator value"""
    b0::Any
    """Index of the fluctuator to be flipped next"""
    next_idx::Any
    """Time interval for next flip event"""
    next_τ::Any
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
struct StochasticNoise{has_sym} <: AbstractOpenSys
    """System-bath coupling operator"""
    ops::AbstractCouplings
    """Symbol used for DEDataArray type"""
    sym::Union{Nothing,Symbol}
end

StochasticNoise(ops::AbstractCouplings, ::Nothing) =
    StochasticNoise{false}(ops, nothing)
StochasticNoise(ops::AbstractCouplings, sym::Symbol) =
    StochasticNoise{true}(ops, sym)

#TODO: OpenSys interface update_cache!, update_vectorized_cache!, QTBase.update_ρ!
(S::StochasticNoise{false})(A, n, tf::Real, t) =
    A .+= -1.0im * tf * sum(n .* S.ops(t))
(S::StochasticNoise{false})(A, n, tf::UnitTime, t) = S(A, n, 1.0, t / tf)
(S::StochasticNoise{true})(A, u::DEDataArray, tf::Real, t) =
    A .+= -1.0im * tf * sum(getfield(u, S.sym) .* S.ops(t))
(S::StochasticNoise{true})(A, u::DEDataArray, tf::UnitTime, t) =
    S(A, u, 1.0, t / tf)

# ============ StochasticNoise{true} =================
QTBase.update_ρ!(
    du,
    u::DEDataArray,
    p::ODEParams,
    t,
    S::StochasticNoise{true},
) = QTBase.update_ρ!(du, u, p.tf, t, S)

function QTBase.update_ρ!(
    du,
    u::DEDataArray,
    tf::Real,
    t,
    S::StochasticNoise{true},
)
    H = sum(getfield(u, S.sym) .* S.ops(t))
    BLAS.gemm!('N', 'N', -1.0im * tf, H, u.x, 1.0 + 0.0im, du.x)
    BLAS.gemm!('N', 'N', 1.0im * tf, u.x, H, 1.0 + 0.0im, du.x)
end

QTBase.update_ρ!(
    du,
    u::DEDataArray,
    tf::UnitTime,
    t,
    S::StochasticNoise{true},
) = QTBase.update_ρ!(du, u, 1.0, t / tf, S)

# ============ StochasticNoise{false} =================
QTBase.update_ρ!(du, u, p::ODEParams, t, S::StochasticNoise{false}) =
    QTBase.update_ρ!(du, u, p.control, p.tf, t, S)
QTBase.update_ρ!(
    du,
    u,
    ctr::FluctuatorControl,
    tf,
    t,
    S::StochasticNoise{false},
) = QTBase.update_ρ!(du, u, ctr(), tf, t, S.ops)

function QTBase.update_ρ!(du, u, n, tf::Real, t, couplings::AbstractCouplings)
    H = sum(n .* couplings(t))
    BLAS.gemm!('N', 'N', -1.0im * tf, H, u, 1.0 + 0.0im, du)
    BLAS.gemm!('N', 'N', 1.0im * tf, u, H, 1.0 + 0.0im, du)
end

function QTBase.update_ρ!(
    du,
    u::DEDataArray,
    n,
    tf::Real,
    t,
    couplings::AbstractCouplings,
)
    H = sum(n .* couplings(t))
    BLAS.gemm!('N', 'N', -1.0im * tf, H, u.x, 1.0 + 0.0im, du.x)
    BLAS.gemm!('N', 'N', 1.0im * tf, u.x, H, 1.0 + 0.0im, du.x)
end

QTBase.update_ρ!(du, n, tf::UnitTime, t, S::StochasticNoise) =
    QTBase.update_ρ!(du, n, 1.0, t / tf, S)

# ============== StochasticNoise{true} ================
QTBase.update_vectorized_cache!(
    cache,
    u::DEDataArray,
    p::ODEParams,
    t,
    S::StochasticNoise{true},
) = QTBase.update_vectorized_cache!(cache, u, p.tf, t, S)

function QTBase.update_vectorized_cache!(
    cache,
    u::DEDataArray,
    tf::Real,
    t,
    S::StochasticNoise{true},
)
    H = sum(getfield(u, S.sym) .* S.ops(t))
    iden = Matrix{eltype(H)}(I, size(H))
    cache .= 1.0im * tf * (transpose(H) ⊗ iden - iden ⊗ H)
end

QTBase.update_vectorized_cache!(
    cache,
    u::DEDataArray,
    tf::UnitTime,
    t,
    S::StochasticNoise{true},
) = QTBase.update_vectorized_cache!(cache, u, 1.0, t / tf, S)

# ============ StochasticNoise{false} =================
QTBase.update_vectorized_cache!(
    cache,
    u,
    p::ODEParams,
    t,
    S::StochasticNoise{false},
) = QTBase.update_vectorized_cache!(cache, p.control, p.tf, t, S)

QTBase.update_vectorized_cache!(
    cache,
    ctr::FluctuatorControl,
    tf,
    t,
    S::StochasticNoise{false},
) = QTBase.update_vectorized_cache!(cache, ctr(), tf, t, S)

function QTBase.update_vectorized_cache!(
    cache,
    n::AbstractArray,
    tf::Real,
    t,
    S::StochasticNoise{false},
)
    H = sum(n .* S.ops(t))
    iden = Matrix{eltype(H)}(I, size(H))
    cache .= 1.0im * tf * (transpose(H) ⊗ iden - iden ⊗ H)
end

QTBase.update_vectorized_cache!(cache, n, tf::UnitTime, t, S::StochasticNoise) =
    QTBase.update_vectorized_cache!(cache, n, S, 1.0, t / tf)
