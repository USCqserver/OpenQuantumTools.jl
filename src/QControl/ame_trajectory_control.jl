"""
$(TYPEDEF)

Defines a adiabatic master equation trajectory controller

# Fields

$(FIELDS)
"""
mutable struct AMETrajectoryControl <: AbstractAnnealingControl
    op
    r::Float64
end
reset!(c::AMETrajectoryControl, u0, initializer) = c.r = rand()


mutable struct AMETrajectoryDEControl <: AbstractAnnealingControl
    op
    sym::Symbol
end


get_controller_name(::Union{AMETrajectoryControl,AMETrajectoryDEControl}) =
    :ame_trajectory_control
reset!(c::AMETrajectoryDEControl, u0, initializer) =
    setfield!(u0, c.sym, rand())

AMETrajectoryControl(op, ::Nothing) = AMETrajectoryControl(op, rand())
AMETrajectoryControl(op, sym::Symbol) = AMETrajectoryDEControl(op, sym)


function build_callback(::AMETrajectoryControl)
    condition = function (u, t, integrator)
        real(u' * u) - integrator.p.control.r
    end

    affect! = function (integrator)
        A = QTBase.ame_jump(
            integrator.p.control.op,
            integrator.u,
            integrator.p.tf,
            integrator.t,
        )
        integrator.p.control.r = rand()
        new_state = normalize(A * integrator.u)
        for c in full_cache(integrator)
            c .= new_state
        end
    end

    ContinuousCallback(condition, affect!)
end


function build_callback(::AMETrajectoryControl, controller_name::Symbol)
    condition = function (u, t, integrator)
        controller = getproperty(integrator.p.control, controller_name)
        real(u' * u) - controller.r
    end

    affect! = function (integrator)
        controller = getproperty(integrator.p.control, controller_name)
        A = QTBase.ame_jump(
            controller.op,
            integrator.u,
            integrator.p.tf,
            integrator.t,
        )
        controller.r = rand()
        new_state = normalize(A * integrator.u)
        for c in full_cache(integrator)
            c .= new_state
        end
    end

    ContinuousCallback(condition, affect!)
end


function build_callback(ctr::AMETrajectoryDEControl)
    condition = function (u, t, integrator)
        real(u.x' * u.x) - getfield(u, integrator.p.control.sym)
    end

    affect! = function (integrator)
        A = QTBase.ame_jump(
            integrator.p.control.op,
            integrator.u,
            integrator.p.tf,
            integrator.t,
        )
        new_r = rand()
        new_state = normalize(A * integrator.u.x)
        for c in full_cache(integrator)
            c.x .= new_state
            setfield!(c, integrator.p.control.sym, new_r)
        end
    end

    ContinuousCallback(condition, affect!)
end


function build_callback(ctr::AMETrajectoryDEControl, controller_name::Symbol)
    condition = function (u, t, integrator)
        controller = getproperty(integrator.p.control, controller_name)
        real(u.x' * u.x) - getfield(u, controller.sym)
    end

    affect! = function (integrator)
        A = QTBase.ame_jump(
            controller.op,
            integrator.u,
            integrator.p.tf,
            integrator.t,
        )
        new_r = rand()
        new_state = normalize(A * integrator.u.x)
        for c in full_cache(integrator)
            c.x .= new_state
            setfield!(c, controller.sym, new_r)
        end
    end

    ContinuousCallback(condition, affect!)
end
